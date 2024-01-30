/*
This is my first attempt at a naive pathtracer. This shader uses Global Illumination
techniques and equations to mimic realistics light behaviour in a given 3d scene.

If you would like to learn more about light rendering and pathtracers I cannot
recommend enough the TU Wien course available for free on Youtube:
https://www.youtube.com/watch?v=pjc1QAI6zS0&list=PLujxSBD-JXgnGmsn7gEyN28P1DnRZG7qi&ab_channel=TwoMinutePapers
This is a big shout-out to TU Wien University and professors Károly Zsolnai and Thomas
Auzinger for making this wonderful course available to everyone on the internet.

This program is heavily based on the basic pathtracer named SmallPaint presented
in the course and proposes a much faster convergence using gpu parallelization
instead of simple cpu multithreading. To learn more you can find its source code here:
https://users.cg.tuwien.ac.at/zsolnai/gfx/smallpaint/

This is a multipass shader that stores for each pixel the sum of all the previous
ray colors that were shot from the pixel. In the alpha channel we store the total
number of samples that were shot so that we can calculate the average of all the
rays color by the number of rays. We do this in the other shader that is the entry
point of the application.

Feel free to reuse any part of the code that you see. Enjoy!

William Cantin
*/

//Here are various uniforms provided by Shadertoy that could be useful
//iResolution, iGlobalTime (also as iTime),
//iTimeDelta, iFrame, iMouse, iMouseButton, iDate, iSampleRate,
//iChannelN with N in [0, 9] and iChannelResolution[] are available uniforms.

#ifdef GL_ES
precision mediump float;
#endif

//#iChannel0"self"
//#iChannel0::MagFilter"Nearest"

//-------------------------------------------------------------------
// Constants
//-------------------------------------------------------------------
const int SAMPLE_PER_PIXEL=1;
const int DEPTH_MAX=10;

const float RUSSIAN_ROULETTE_FACTOR=1.1111111111;
const float RUSSIAN_ROULETTE_STOP_PROBABILITY=.1;

const int LIGHT=0;
const int BRDF_DIFFUSE=1;
const int BRDF_SPEC=2;
const int BRDF_REFR=3;

const float REFR_INDEX_GLASS=1.5;
const float REFR_INDEX_AIR=1.;

const float PI=3.1415926536;
const float EPS=1e-6;
const float INF=1e9;

//-------------------------------------------------------------------
// Structs
//-------------------------------------------------------------------
struct Camera{
    vec3 pos,lookAt;
};

struct Material{
    vec3 color;
    float emission;
    int typeBrdf;
};

struct Sphere{
    vec3 center;
    float radius;
    Material mat;
};

struct Plane{
    vec3 normal;
    float d;
    Material mat;
};

struct Ray{
    vec3 origin;
    vec3 direction;
};

struct Intersection{
    float closest;
    Material mat;
    vec3 pos;
    vec3 normal;
};

//-------------------------------------------------------------------
// Helper functions
//-------------------------------------------------------------------

//This is a random number generator between 0 and 1
// All credits to Inigo Quilez that provided this function
float seed;//Defined in the main function
float randomGenerator0(){
    return fract(sin(seed++)*43758.5453123);
}

// This is a random number generator between -1 and 1
float randomGenerator1(){
    return((randomGenerator0()*2.)-1.);
}

//Generate a random point located on the hemisphere of radius 1
// located at (0,0,0)
vec3 hemisphere()
{
    float u1=randomGenerator0();
    
    float r=sqrt(1.-u1*u1);
    float phi=2.*PI*randomGenerator0();
    return vec3(cos(phi)*r,sin(phi)*r,u1);
}

//Create an orthonormal basis from only 1 vector named v1
// (we assume v1 is already normalized)
mat3 createOrthormalBase(vec3 v1){
    mat3 mat=mat3(0.);
    
    //If we are completely aligned with the y axis it is trivial to find the
    // other vectors (its z and x)
    if(v1.y==1.||v1.y==-1.){
        mat=mat3(v1,vec3(1.,0.,0.),vec3(0.,0.,1.));
    }
    
    else{
        //We find an orthogonal vector with v1 and a random up direction vector
        vec3 v2=normalize(cross(v1,vec3(0.,1.,0.)));
        
        //We find another vector orthogonal to the first 2
        vec3 v3=normalize(cross(v1,v2));
        mat=mat3(v1,v2,v3);
    }
    
    return mat;
}

//-------------------------------------------------------------------
// Plane objects
//-------------------------------------------------------------------
float intersectionPlane(Plane plane,Ray ray){
    //D0 is the cos between ray direction and normal of the plane
    float d0=dot(plane.normal,ray.direction);
    
    //If the ray is not parallel to the plane
    if(d0!=0.)
    {
        float t=-1.*(((dot(plane.normal,ray.origin))+plane.d)/d0);
        
        //If t < eps, it would mean that the ray is already starting from the plane
        // so we ignore it
        return(t>EPS)?t:0.;
    }
    
    // If we are parallel to the plane return 0 (we dont hit the plane)
    else{
        return 0.;
    }
}

//-------------------------------------------------------------------
// Sphere object
//-------------------------------------------------------------------

//Return the normal of the sphere from a certain point on the surface
vec3 getSphereNormal(Sphere sphere,vec3 pos){
    return normalize(pos-sphere.center);
}

//Return the closest intersection of the sphere and the origin of the ray. We ignore
// rays that start from the sphere (distance 0) and the intersections behind the ray (negative number)
float intersectionSphere(Sphere sphere,Ray ray){
    //To find the intersection, we use a polynomial function and use the paramtric equation of a sphere
    // For more info you can learn about it in more detail in the Unit 02 - The Rendering Equation and
    // Ray Tracing (Károly Zsolnai-Fehér) powerpoint presentation at this page
    // https://www.cg.tuwien.ac.at/courses/Rendering/VU.SS2019.html
    
    float b=dot(((ray.origin-sphere.center)*2.),(ray.direction));
    float c=dot((ray.origin-sphere.center),
    (ray.origin-sphere.center))
    -(sphere.radius*sphere.radius);
    
    //Disc is the discrimant in the polynomial function (-b +/- sqrt(b*b - 4*a*c ) ) 2a
    // We dont calculate the square root yet for optimisation purposes
    float discriminant=(b*b)-(4.*c);
    
    //If discrimant is negative there are no intersections
    // (sqrt of negative number is impossible here)
    if(discriminant<0.){
        return 0.;
    }
    else{
        discriminant=sqrt(discriminant);
    }
    
    //We dont yet divide by 2 for optimisation reasons
    float sol1=-b+discriminant;
    float sol2=-b-discriminant;
    
    //We return the closest solution that is over epsilon.
    // If we were < epsilon, it would mean
    // the ray is originated from the sphere so we ignore it
    return(sol2>EPS)?(sol2/2.):((sol1>EPS)?sol1/2.:0.);
}

//-------------------------------------------------------------------
// Scene creation
//-------------------------------------------------------------------
const Sphere lights[]=Sphere[](
    Sphere(vec3(0.,-5.,-4.),2.,Material(vec3(1.,1.,1.),10.,LIGHT))
);

//center, radius, Mat (color, emission, type_brdf)
const Sphere worldSphere[]=Sphere[](
    Sphere(vec3(0.,0,-4.),1.,Material(vec3(0.,1.,1.),0.,BRDF_DIFFUSE)), // Teal ball
    Sphere(vec3(0.,-1.5,-4.),.5,Material(vec3(1.,.8863,.7373),0.,BRDF_REFR)),
    Sphere(vec3(1.,-1.5,-4.),.25,Material(vec3(1.,.8863,.7373),0.,BRDF_SPEC)),
    Sphere(vec3(-2.,0,-4.),1.,Material(vec3(1.,.8863,.7373),0.,BRDF_SPEC)),
    Sphere(vec3(2.,0,-4.),1.,Material(vec3(1.,.8863,.7373),0.,BRDF_REFR)),
    Sphere(vec3(0.,0.,-1.5),.5,Material(vec3(1.,.8863,.7373),0.,BRDF_REFR))
);

//normal, d
const Plane worldPlane[]=Plane[](
    Plane(vec3(0.,-1.,0.),1.,Material(vec3(1.,1.,1.),0.,BRDF_DIFFUSE)),//floor
    Plane(vec3(0.,0.,1),5.5,Material(vec3(1.,1.,1.),0.,BRDF_DIFFUSE)),// front of cam
    Plane(vec3(-1.,0.,0.),5.5,Material(vec3(.2824,1.,0.),0.,BRDF_DIFFUSE)),//right
    Plane(vec3(1.,0.,0.),5.5,Material(vec3(1.,0.,0.),0.,BRDF_DIFFUSE)),// left
    //Plane(vec3(0.,1.,0.),6.,Material(vec3(1.,1.,1.),0.,BRDF_DIFFUSE)),//Roof
    Plane(vec3(0.,0.,-1.),5.,Material(vec3(.5216,.5216,.5216),0.,BRDF_DIFFUSE))//Back of cam
);

//Find the nearest object hit by the ray
Intersection intersectionScene(Ray ray){
    float closest=INF;
    vec3 hitPosition=vec3(INF,INF,INF);
    Material mat=Material(vec3(0.),0.,BRDF_DIFFUSE);
    vec3 normal=vec3(0.,0.,0.);
    
    //Spheres intersection
    for(int i=0;i<worldSphere.length();i++){
        float distance=intersectionSphere(worldSphere[i],ray);
        if(distance<closest&&distance>EPS){
            closest=distance;
            mat=worldSphere[i].mat;
            hitPosition=ray.origin+ray.direction*closest;
            normal=getSphereNormal(worldSphere[i],hitPosition);
        }
    }
    
    //Planes
    for(int i=0;i<worldPlane.length();i++){
        float distance=intersectionPlane(worldPlane[i],ray);
        if(distance<closest&&distance>EPS){
            closest=distance;
            mat=worldPlane[i].mat;
            hitPosition=ray.origin+ray.direction*closest;
            normal=worldPlane[i].normal;
        }
    }
    
    //Light
    for(int i=0;i<lights.length();i++){
        float distance=intersectionSphere(lights[i],ray);
        if(distance<closest&&distance>EPS){
            closest=distance;
            mat=lights[i].mat;
            hitPosition=ray.origin+ray.direction*closest;
            normal=getSphereNormal(lights[i],hitPosition);
        }
    }
    
    return Intersection(closest,mat,hitPosition,normal);
}

//-------------------------------------------------------------------
// Camera
//-------------------------------------------------------------------

/*
Big thanks to the youtube channel "The Art of Code" for showing me this method.
You can learn more about cameras here:
https://www.youtube.com/watch?v=PBxuVlp7nuM&ab_channel=TheArtofCode
*/
Ray createRayFromCamera(vec2 uv){
    // Normalized pixel coordinates (from 0 to 1)
    uv-=.5;
    uv.x*=iResolution.x/iResolution.y;
    
    float w=iResolution.x;
    float h=iResolution.y;
    
    //We add some randomness to each pixel so that we dont always shoot
    // through its center
    uv.x=uv.x+randomGenerator1()/w;
    uv.y=uv.y+randomGenerator1()/h;
    
    //We create a camera with a position and its lookat position
    Camera camera=Camera(vec3(0.,-0.,4.),vec3(0.));
    
    vec3 forward=normalize(camera.lookAt-camera.pos);
    vec3 right=normalize(cross(forward,vec3(0.,1.,0.)));
    
    //No need to normalize here because
    //forward and right already normalized
    vec3 up=(cross(forward,right));
    
    float zoom=1.;
    vec3 centerOfPlane=camera.pos+forward*zoom;
    
    Ray ray=Ray(camera.pos,normalize((centerOfPlane+right*uv.x+up*uv.y)-camera.pos));
    
    return ray;
}

//-------------------------------------------------------------------
// Path tracing
//-------------------------------------------------------------------
vec3 trace(Ray ray){
    bool continueTracing=true;
    vec3 color=vec3(0);
    vec3 factor=vec3(1.);
    int depth=0;
    
    while(continueTracing){
        
        // Russian roulette: starting at depth 5, each recursive step will stop with a certain probability
        
        if(depth>DEPTH_MAX)
        {
            
            //We use Russian roulette to terminate path and cull
            // the least interesting rays. You can learn more about it
            // in the wonderful explanation written by user RichieSams on
            // stack exchange: https://computergraphics.stackexchange.com/questions/2316/is-russian-roulette-really-the-answer
            
            // //We stop according to the probability
            // if(randomGenerator0()<=RUSSIAN_ROULETTE_STOP_PROBABILITY)
            // {
                //     continueTracing=false;
                //     return color;
            // }
            
            continueTracing=false;
            return vec3(0.);
        }
        else{
            //We find the closest intersection with the ray
            Intersection intersection=intersectionScene(ray);
            
            //If we dont hit anything we are done
            if(intersection.closest==INF){
                continueTracing=false;
                return vec3(0.);
            }
            
            // Travel the ray to the hit point where the closest object lies and
            // compute the surface normal there.
            
            //We set the new ray's origin to the coordinates we found
            // We add a small epsilon to the position to make sure we dont run
            // into Surface Acne. Thanks to user Jon Hanson on stack overflow:
            // https://stackoverflow.com/questions/47949652/how-can-i-improve-my-path-tracer
            ray.origin=intersection.pos+(.0001)*intersection.normal;
            
            if(intersection.mat.typeBrdf==LIGHT){
                // color=intersection.mat.color*intersection.mat.emission;
                color=factor*intersection.mat.color*intersection.mat.emission;
                return color;
            }
            
            else if(intersection.mat.typeBrdf==BRDF_DIFFUSE){
                //We pick a random direction on the hemisphere of reflection
                vec3 sampledDir=hemisphere();
                
                //We create an orthonormal basis with the normal of the surface hit
                mat3 orthonormalBasis=createOrthormalBase(intersection.normal);
                vec3 rotX=orthonormalBasis[1];
                vec3 rotY=orthonormalBasis[2];
                
                //We calculate the new ray Direction according to the orthonormal base
                vec3 rotatedDir;
                rotatedDir.x=dot(vec3(rotX.x,rotY.x,intersection.normal.x),sampledDir);
                rotatedDir.y=dot(vec3(rotX.y,rotY.y,intersection.normal.y),sampledDir);
                rotatedDir.z=dot(vec3(rotX.z,rotY.z,intersection.normal.z),sampledDir);
                
                ray.direction=rotatedDir;//already normalized
                
                float cosTheta=dot(ray.direction,intersection.normal);
                
                factor=intersection.mat.color*cosTheta*factor;
                
                color=intersection.mat.color*factor*color;
                // //__________________________________________________
                // vec3 lightVector=normalize(lights[0].center-intersection.pos);
                // float cosTheta2=dot(intersection.normal,lightVector);
                
                // color=cosTheta2*intersection.mat.color;
                
                // return color;
                // //__________________________________________________
                
            }
            
            else if(intersection.mat.typeBrdf==BRDF_SPEC){
                ray.direction=reflect(ray.direction,intersection.normal);
            }
            
            // Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
            // to compute the outgoing reflection and refraction directions and probability weights.
            // For more info on the vector form of Snell's law you can visit this page:
            // https://en.wikipedia.org/wiki/Snell%27s_law
            else if(intersection.mat.typeBrdf==BRDF_REFR){
                float n1=REFR_INDEX_AIR;
                float n2=REFR_INDEX_GLASS;
                
                //Note the minus behind the normal to fit Sclick's equation
                float cosTheta1=dot(-intersection.normal,ray.direction);
                
                //If the direction of the ray and the normal are in the same direction,
                // we are inside the glass medium so we need to revert the normal
                // and swap the refractive indexes. We dont need to recalculate the new
                // cosTheta and can just make it positive.
                if(cosTheta1<0.){
                    n1=REFR_INDEX_GLASS;
                    n2=REFR_INDEX_AIR;
                    intersection.normal=-intersection.normal;
                    cosTheta1=-cosTheta1;
                    
                }
                
                //R0 is the probability of reflection at angle 0 (perfectly aligned with the normal)
                // It is going to be useful when using Schlicks approximation
                float R0=(n1-n2)/(n1+n2);
                R0=R0*R0;
                
                //We calculate the probability of reflection for the angle we have
                float reflectionProbability=R0+(1.-R0)*pow(1.-cosTheta1,5.);
                vec3 reflectedDirection=reflect(ray.direction,intersection.normal);
                
                float ratioRefractive=n1/n2;
                
                //We find cosTheta2. We dont do the square root yet for optimisation reasons
                // in case the term is negative (sqrt of negative is not possible here)
                float cosTheta2=1.-ratioRefractive*ratioRefractive*(1.-cosTheta1*cosTheta1);
                
                //We make sure cosTheta2 is positive for the sqrt and we randomly send the ray
                // according to the probability of reflection
                if(cosTheta2>0.&&randomGenerator0()>reflectionProbability){
                    ray.direction=normalize(ratioRefractive*ray.direction+
                        (ratioRefractive*cosTheta1-sqrt(cosTheta2))*intersection.normal);
                        
                        //Special case here where the position of the new ray
                        // depends on which side of the surface we are going to
                        
                        //If we stay inside the medium
                        if(dot(ray.direction,intersection.normal)>0.){
                            ray.origin=intersection.pos+(.0001)*intersection.normal;
                        }
                        
                        //We leave the medium
                        else{
                            ray.origin=intersection.pos-(.0001)*intersection.normal;
                        }
                    }
                    
                    //We reflect the ray
                    else{
                        ray.direction=reflectedDirection;
                    }
                }
                
            }
            
            depth++;
        }
        
        return color;
    }
    
    //-------------------------------------------------------------------
    // Main function
    //-------------------------------------------------------------------
    void mainImage(out vec4 fragColor,in vec2 fragCoord)
    {
        //We make sure our coordinates are between 0 and 1
        vec2 uv=fragCoord/iResolution.xy;
        
        //Seed our random number generator with dynamic value like itime
        seed=iTime+iResolution.y*fragCoord.x/iResolution.x+fragCoord.y/iResolution.y;
        
        vec3 color=vec3(0.);
        
        //We trace SAMPLE_PER_PIXEL rays for each pixel and sum their returned colors
        for(int i=0;i<SAMPLE_PER_PIXEL;i++){
            //Ray ray=Ray(vec3(0,0,0),normalize(cam-vec3(0,0,0)));
            Ray ray=createRayFromCamera(uv);
            
            // Ray ray=Ray(vec3(fragCoord,-1.),vec3(0.,0.,1.));
            Intersection intersection=intersectionScene(ray);
            
            color+=trace(ray);
        }
        
        //We clamp the color to make sure it is positive (in case of
        // float precision errors)
        color=clamp(color,0.,INF);
        
        vec4 lastFrameColor=texture(iChannel0,uv).rgba;
        
        //The output texture of this shader stores the sum of the
        // radiance of all the rays traced per pixel and the total number
        // of samples taken in its alpha channel
        fragColor=lastFrameColor+vec4(color,SAMPLE_PER_PIXEL);
    }
    
