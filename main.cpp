#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "Vector3f.h"
#include <math.h>

using namespace parser;
using namespace std;
typedef unsigned char RGB[3];

Vector3f k;
double tdistance;
int mid;
int VvsS = 0;
int point = 0;
int Meshface = 0;
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
    if(delta < shadow_ray_epsilon)
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
    
    ray.cam.x = camera.position.x;
    ray.cam.y = camera.position.y;
    ray.cam.z = camera.position.z;

    ray.direction.x = u;
    ray.direction.y = v;
    Vector3f a (camera.gaze.x,camera.gaze.y,camera.gaze.z);
    Vector3f NorGaze = a.normalize();
    Vector3f v1(u+camera.near_distance*NorGaze.x,-v-camera.near_distance*NorGaze.y,camera.near_distance*NorGaze.z);
    ray.direction = v1;
    return ray;
}
Vector3f make_color( Scene scene,Ray ray ){

    
    Vector3f AL(scene.ambient_light.x,scene.ambient_light.y,scene.ambient_light.z);                                             // ambient_light
    Vector3f D(scene.materials[mid].diffuse.x,scene.materials[mid].diffuse.y,scene.materials[mid].diffuse.z);                   // diffuse coficient
    Vector3f PLP(scene.point_lights[0].position.x,scene.point_lights[0].position.y,scene.point_lights[0].position.z);                 // light position
    Vector3f MA(scene.materials[mid].ambient.x,scene.materials[mid].ambient.y,scene.materials[mid].ambient.z);                  // material ambiant
    Vector3f PLI(scene.point_lights[0].intensity.x,scene.point_lights[0].intensity.y,scene.point_lights[0].intensity.z);              // ligth instensity     
    Vector3f S(scene.materials[mid].specular.x,scene.materials[mid].specular.y,scene.materials[mid].specular.z);                // specular
    Vector3f DirecDistan(tdistance*ray.direction.x ,tdistance*ray.direction.y,tdistance*ray.direction.z);
    Vector3f s = (ray.cam + DirecDistan);
    Vector3f n;
    if (VvsS == 2 ){ // sphere
        Vec3f test_r = scene.vertex_data[scene.spheres[point].center_vertex_id-1];
        Vector3f c_sphere(test_r.x,test_r.y,test_r.z);
        n.x = (s.x-c_sphere.x)/scene.spheres[point].radius;
        n.y = (s.y-c_sphere.y)/scene.spheres[point].radius;
        n.z = (s.z-c_sphere.z)/scene.spheres[point].radius;
        //n = n.normalize();
        cout<< tdistance<<" " << VvsS << endl;
    }
    if (VvsS == 1 ){ // triangle

        Vector3f a (scene.vertex_data[scene.triangles[point].indices.v0_id-1].x,scene.vertex_data[scene.triangles[point].indices.v0_id-1].y,scene.vertex_data[scene.triangles[point].indices.v0_id-1].z);
        Vector3f b (scene.vertex_data[scene.triangles[point].indices.v1_id-1].x,scene.vertex_data[scene.triangles[point].indices.v1_id-1].y,scene.vertex_data[scene.triangles[point].indices.v1_id-1].z);
        Vector3f c (scene.vertex_data[scene.triangles[point].indices.v2_id-1].x,scene.vertex_data[scene.triangles[point].indices.v2_id-1].y,scene.vertex_data[scene.triangles[point].indices.v2_id-1].z);
        n = ((b - a).cross(c - a )).normalize(); 

    }
    if (VvsS == 3){ // mesh

        Vector3f a (scene.vertex_data[scene.meshes[point].faces[Meshface].v0_id-1].x,scene.vertex_data[scene.meshes[point].faces[Meshface].v0_id-1].y,scene.vertex_data[scene.meshes[point].faces[Meshface].v0_id-1].z);
        Vector3f b (scene.vertex_data[scene.meshes[point].faces[Meshface].v1_id-1].x,scene.vertex_data[scene.meshes[point].faces[Meshface].v1_id-1].y,scene.vertex_data[scene.meshes[point].faces[Meshface].v1_id-1].z);
        Vector3f c (scene.vertex_data[scene.meshes[point].faces[Meshface].v2_id-1].x,scene.vertex_data[scene.meshes[point].faces[Meshface].v2_id-1].y,scene.vertex_data[scene.meshes[point].faces[Meshface].v2_id-1].z);
        n = ((b - a).cross(c - a )).normalize(); 
    }
    Vector3f l = ( PLP - s ).normalize();
    Vector3f v = ( ray.cam - s ).normalize();
    Vector3f h = (l + v ).normalize();
    double distance = (PLP - s ).length();
    double nl = n.dot(l);
    double hn = h.dot(n);
    double maxnl = 0;
    double maxhn = 0;
    if (nl > 0)
        maxnl =10* nl;
    if ( hn > 0 )
        maxhn =10* hn;
    Vector3f L = AL * MA + ( D * PLI * maxnl /  distance * distance ) +  ( S * PLI * pow(maxhn,scene.materials[mid].phong_exponent) / distance * distance );

    return L;
    
}
float clamp(float x, float a, float b){

    return x < a ? a : (x > b ? b : x);
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

    int width = scene.cameras[0].image_width, height = scene.cameras[0].image_height;
    int  i = 0;
    Camera camera;
    Ray ray;
    unsigned char* image = new unsigned char [width * height * 3];
    
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {   
            ray = raytracer(x,y,width,height,scene.cameras[0]);
            int j = 0 ;
            while ( j < scene.spheres.size()){
                if (  sphere ( ray,scene.spheres[j],scene.vertex_data,scene.shadow_ray_epsilon) ){
                    mid = scene.spheres[j].material_id-1;
                    VvsS = 2 ;
                    point = j;
                    break;
                }
                j++;
            }
            j = 0 ;

            while ( j < scene.triangles.size()){
                
                if( triangle(ray, scene.triangles[j], scene.vertex_data) ){
                    mid = scene.triangles[j].material_id-1;
                    VvsS = 1 ;
                    point = j;
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
                        mid = triangle_temp.material_id-1;
                        VvsS = 3 ;
                        point = j;
                        Meshface = k;
                        break;
                    }
                    k++;
                }
                if(VvsS == 3)
                        break;
                j++;
            } 
            
            if (tdistance > 0){

                Vector3f color = make_color( scene,ray);

                image[i++] = clamp(int(color.x), 0, 255);
                image[i++] = clamp(int(color.y), 0, 255);
                image[i++] = clamp(int(color.z), 0, 255);
                //cout << int(color.x) % 256  << " " << int(color.y)  % 256<< " " << int(color.z)  % 256<< endl;
                

            }else{
                image[i++] = scene.background_color.x;
                image[i++] = scene.background_color.y;
                image[i++] = scene.background_color.z;

            }
            tdistance = 0;
            VvsS = 0 ;
            point = 0;
            Meshface = 0;
        }
    }
    //cout<< j << endl;
    write_ppm(argv[2], image, width, height);

}