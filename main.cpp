#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "Vector3f.h"
using namespace parser;
using namespace std;
typedef unsigned char RGB[3];

struct Ray{

    Vector3f cam;
    Vector3f direction;
};

bool sphere (Ray ray ,Sphere sphere , std::vector<Vec3f> vertex_data ){
    Vec3f test_r = vertex_data[sphere.center_vertex_id-1];
    Vector3f c_sphere(test_r.x,test_r.y,test_r.z);
    ray.direction.z = test_r.z;
    //ray.direction.normalize();
    //Vector3f d = ray.direction - c_sphere ;
    //Vector3f ec = ray.cam - c_sphere ;
    /*float bb = 2 * d.dot(ec) ;
    float a = d.dot(d);
    float c = ec.dot(ec);
    if ( (bb*bb - 4 * a*c ) >= 0 )
        return false;
    return true;*/
    double a,b,c,delta,t,t1,t2;
    double c1=(ray.cam.x-c_sphere.x)*(ray.cam.x-c_sphere.x)+(ray.cam.y-c_sphere.y)*(ray.cam.y-c_sphere.y)+(ray.cam.z-c_sphere.z)*(ray.cam.z-c_sphere.z);//.dot(ray.cam-c_sphere);
    c=c1-sphere.radius*sphere.radius;
    b=2*(ray.direction.x)*(ray.cam.x-c_sphere.x)+2*(ray.direction.y)*(ray.cam.y-c_sphere.y)+2*(ray.direction.z)*(ray.cam.z-c_sphere.z);
    a=(ray.direction.x*ray.direction.x)+(ray.direction.y*ray.direction.y)+(ray.direction.z*ray.direction.z);
    delta=b*b-4*a*c;
    /*Vector3f r1 = c_sphere-ray.direction;
    float lr1 = r1.length();
    std:: cout << lr1 << endl;
    if (lr1 <= sphere.radius)
        return true;
    return false;*/
    if(delta<1e-3)
    return false;
    else
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
    float M =  ( a1 * ( e * i - h * f ) + b1 * ( g * f - d1 * i ) + c1 * ( d1 * h - g * e) );
    float beta =  ( j * ( e * i - h * f ) + k * ( g * f - d1 * i ) + l * ( d1 * h - g * e) ) / M;
    float alfa =  ( i * ( a1 * k - j * b1 ) + h * ( j * c1 - a1 * l ) + g * ( b1 * l - k * c1 ) ) / M;
    float t    = - ( f * ( a1 * k - j * b1 ) + e * ( j * c1 - a1 * l ) + d1 * ( b1 * l - k * c1 ) ) / M;
    if ( alfa < 0  || alfa > 1)
        return false;
    if (beta < 0 || beta > 1 - alfa )
        return false;
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
    //ray.direction.z = camera.near_distance;
    Vector3f v1(u+camera.near_distance*camera.gaze.x,v+camera.near_distance*camera.gaze.y,-camera.near_distance*camera.gaze.z);
    //ray.direction=u+v;
    ray.direction = v1;
    //ray.direction.normalize();
    return ray;
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

    int width = 800, height = 800;
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
            if (  sphere ( ray,scene.spheres[0],scene.vertex_data) || triangle(ray, scene.triangles[0], scene.vertex_data) ){
                image[i++] = BAR_COLOR[0][0];
                image[i++] = BAR_COLOR[0][1];
                image[i++] = BAR_COLOR[0][2];
            }else{
                image[i++] = BAR_COLOR[7][0];
                image[i++] = BAR_COLOR[7][1];
                image[i++] = BAR_COLOR[7][2];
            }

        }
    }

    write_ppm(argv[2], image, width, height);

}
