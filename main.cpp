#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "Vector3f.h"
#include <math.h>
using namespace parser;
using namespace std;
typedef unsigned char RGB[3];

int tdistance;
int mid;

struct Ray{

    Vector3f cam;
    Vector3f direction;
};

bool sphere (Ray ray ,Sphere sphere , std::vector<Vec3f> vertex_data,float shadow_ray_epsilon ){
    Vec3f test_r = vertex_data[sphere.center_vertex_id-1];
    Vector3f c_sphere(test_r.x,test_r.y,test_r.z);
    double a,b,c,delta;
    double c1=(ray.cam.x-c_sphere.x)*(ray.cam.x-c_sphere.x)+(ray.cam.y-c_sphere.y)*(ray.cam.y-c_sphere.y)+(ray.cam.z-c_sphere.z)*(ray.cam.z-c_sphere.z);//.dot(ray.cam-c_sphere);
    c=c1-sphere.radius*sphere.radius;
    b=2*(ray.direction.x)*(ray.cam.x-c_sphere.x)+2*(ray.direction.y)*(ray.cam.y-c_sphere.y)+2*(ray.direction.z)*(ray.cam.z-c_sphere.z);
    a=(ray.direction.x*ray.direction.x)+(ray.direction.y*ray.direction.y)+(ray.direction.z*ray.direction.z);
    delta=b*b-4*a*c;
    if(delta<shadow_ray_epsilon)
        return false;
    double t1 = -((2*(ray.direction.x)*(ray.cam.x-c_sphere.x)+2*(ray.direction.y)*(ray.cam.y-c_sphere.y)+2*(ray.direction.z)*(ray.cam.z-c_sphere.z))/2 + sqrt( b*b-4*a*c ) )/ 
    ((ray.direction.x*ray.direction.x)+(ray.direction.y*ray.direction.y)+(ray.direction.z*ray.direction.z)) ;
    double t2 = -((2*(ray.direction.x)*(ray.cam.x-c_sphere.x)+2*(ray.direction.y)*(ray.cam.y-c_sphere.y)+2*(ray.direction.z)*(ray.cam.z-c_sphere.z))/2 - sqrt( b*b-4*a*c ) )/ 
    ((ray.direction.x*ray.direction.x)+(ray.direction.y*ray.direction.y)+(ray.direction.z*ray.direction.z)) ;
    
    if ( t1 < t2)
        tdistance = t1;
    else 
        tdistance = t2;

    mid=sphere.material_id;
    
    return true;
}
bool triangle( Ray ray,Triangle triangle, std::vector<Vec3f> vertex_data) {
    Vector3f a (vertex_data[triangle.indices.v0_id-1].x,vertex_data[triangle.indices.v0_id-1].y,vertex_data[triangle.indices.v0_id-1].z);
    Vector3f b (vertex_data[triangle.indices.v1_id-1].x,vertex_data[triangle.indices.v1_id-1].y,vertex_data[triangle.indices.v1_id-1].z);
    Vector3f c (vertex_data[triangle.indices.v2_id-1].x,vertex_data[triangle.indices.v2_id-1].y,vertex_data[triangle.indices.v2_id-1].z);
    Vector3f d = ray.direction - ray.cam;
    float a1 = a.x - b.x;       float b1 = a.y - b.y;       float c1 = a.z - b.z;
    float d1 = a.x - c.x;       float e  = a.y - c.y;       float f  = a.z - c.z;
    float j  = a.x - ray.cam.x; float k  = a.y - ray.cam.y; float l  = a.z - ray.cam.z;
    float g  = d.x;             float h  = d.y;             float i  = d.z;
    float M  =   ( a1 * ( e * i - h * f ) + b1 * ( g * f - d1 * i ) + c1 * ( d1 * h - g * e) );
    float t  = - ( f * ( a1 * k - j * b1 ) + e * ( j * c1 - a1 * l ) + d1 * ( b1 * l - k * c1 ) ) / M;
    double t1 = INFINITY;
    if (tdistance > 0)
        t1 = tdistance;
    if (t < 0 || t > t1)
      return false;
    float alfa =  ( i * ( a1 * k - j * b1 ) + h * ( j * c1 - a1 * l ) + g * ( b1 * l - k * c1 ) ) / M;
    
    if ( alfa < 0  || alfa > 1)
        return false;
    float beta =  ( j * ( e * i - h * f ) + k * ( g * f - d1 * i ) + l * ( d1 * h - g * e) ) / M;
    if (beta < 0 || beta > 1 - alfa )
        return false;
    tdistance=t;
    
    mid=triangle.material_id;
    return true;
}
Ray raytracer(int x, int y , int width, int height, Camera camera){
    Ray ray;
    float u  = (camera.near_plane.x + (camera.near_plane.y-camera.near_plane.x)*(x + 0.5 ) / (width)) ;
    float v  = (camera.near_plane.z + (camera.near_plane.w-camera.near_plane.z)*(y + 0.5 ) / (height));
    float d  =  camera.position.z;
    ray.cam.x = camera.position.x;
    ray.cam.y = camera.position.y;
    ray.cam.z = camera.position.z;

    ray.direction.x = u;
    ray.direction.y = v;
    
    Vector3f v1(u+camera.near_distance*camera.gaze.x,-v-camera.near_distance*camera.gaze.y,camera.near_distance*camera.gaze.z);
    ray.direction = v1;
    return ray;
}
Vector3f make_color(Material material,Vec3f ambient_light,PointLight pointlight,Ray ray){
    Vector3f Vambientlight(ambient_light.x,ambient_light.y,ambient_light.z);                            // ambient_light
    Vector3f Vdiffuse(material.diffuse.x,material.diffuse.y,material.diffuse.z);                        // diffuse coficient
    Vector3f Vpointlightposition(pointlight.position.x,pointlight.position.y,pointlight.position.z);    // light position
    Vector3f VMaterialAmbient(material.ambient.x,material.ambient.y,material.ambient.z);                // material ambiant
    Vector3f VpointlightIntensity(pointlight.intensity.x,pointlight.intensity.y,pointlight.intensity.z);// ligth instensity     
    Vector3f s=ray.direction+ray.cam;                                                                   // s
    Vector3f ls=s-Vpointlightposition;                                                                  // ligth direction
    Vector3f Vspecular(material.specular.x,material.specular.y,material.specular.z);                    // specular

    double dotnl=(s.normalize()).dot(ls.normalize());
    Vector3f h=(ls+( ray.cam - s) ).normalize();
    double dotnh=(s.normalize()).dot(h);
    double max1,max2;
    if(dotnl>0)
        max1=dotnl;
    else
        max1=0;
    if(dotnh>0)
        max2=dotnh;
    else
        max2=0;
   
    Vector3f color=VMaterialAmbient*Vambientlight+Vdiffuse*VpointlightIntensity*max1+Vspecular*VpointlightIntensity*pow(max2,material.phong_exponent);
    return color;
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int width = scene.cameras[0].image_width, height = scene.cameras [0].image_height;
    int columnWidth = width ;
    Camera camera;
    Ray ray;
    unsigned char* image = new unsigned char [width * height * 3];
    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            ray = raytracer(x,y,width,height,scene.cameras[0]);
            int j = 0 ;
            while ( j < scene.spheres.size()){
                if (  sphere ( ray,scene.spheres[j],scene.vertex_data,scene.shadow_ray_epsilon) ){
                    break;
                }
                j++;
            }
            /*j = 0 ;

            while ( j < scene.triangles.size()){
                
                if( triangle(ray, scene.triangles[j], scene.vertex_data) ){
                    break;
                }
                j++;
            }  
                
            j = 0 ;
            while ( j < scene.meshes.size()){
                int k = 0;
                Triangle triangle_temp;
                while(k < scene.meshes[j].faces.size()){
                        
                    triangle_temp.indices = scene.meshes[j].faces[k];
                    triangle_temp.material_id = scene.meshes[j].material_id;
                    if (   triangle(ray, triangle_temp, scene.vertex_data) ){
                        break;
                    }
                    k++;
                }
                j++;
            } 
            */
            if (tdistance > 0){
                Vector3f color = make_color( scene.materials[mid-1],scene.ambient_light,scene.point_lights[0],ray);
               cout << int(color.x) << " " << int(color.y) << " " << int(color.z) << endl;

                image[i++] = int(color.x) % 256;
                image[i++] = int(color.y) % 256;
                image[i++] = int(color.z) % 256;
                //cout << int(color.x) % 256  << " " << int(color.y)  % 256<< " " << int(color.z)  % 256<< endl;
                

            }else{
                image[i++] = scene.background_color.x;
                image[i++] = scene.background_color.y;
                image[i++] = scene.background_color.z;

            }
            tdistance = 0;
        }
    }
    //cout<< j << endl;
    write_ppm(argv[2], image, width, height);

}