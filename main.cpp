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
    ray.direction.normalize();
    Vector3f d = ray.direction - c_sphere ;
    Vector3f ec = ray.cam - c_sphere ;
    float bb = 2 * d.dot(ec) ;
    float a = d.dot(d);
    float c = ec.dot(ec);
    if ( (bb*bb - 4 * a*c ) >= 0 )
        return false;
    return true;
   /* Vector3f r1 = c_sphere-ray.direction;
    float lr1 = r1.length();
    std:: cout << r1.x << endl;
    if (lr1 <= sphere.radius)
        return true;
    return false;*/
}
Ray raytracer(int x, int y , int width, int height, Camera camera){
    Ray ray;
    float u  = -1 + 2*(x + 0.5 ) / width;
    float v  = -1 + 2*(y + 0.5 ) / height;
    float d  =  camera.position.z;
    ray.cam.x = camera.position.x;
    ray.cam.x = camera.position.y;
    ray.cam.x = camera.position.z;
    ray.direction.x = u;
    ray.direction.y = v;
    ray.direction.z = -d;
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

    int width = 640, height = 480;
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
            if ( sphere ( ray,scene.spheres[0],scene.vertex_data)){
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
