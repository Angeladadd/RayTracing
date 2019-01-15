//
//  main.cpp
//  raytracing
//
//  Created by sunchenge on 2018/12/28.
//  Copyright © 2018 sunchenge. All rights reserved.
//

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <sstream>

#include "GLVector.h"
#include "GLObject.h"

#define WINDOW_WIDTH 1200
#define WINDOW_HEIGHT 1200
#define MAX_DEPTH 3
#define AmbientStrength 0.3f
#define AmbientLightColor WHITE


float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

Vec3f trace(Ray &ray, std::vector<Object *> objects, std::vector<Sphere *> lightings, Model *model, int depth){
    float t_min = INFINITY;
    Object *object = NULL;
    Vec3f color = BLACK;
    float bias = 1e-4;
    
    std::size_t obj_size = objects.size();
    std::size_t light_size = lightings.size();
    
    for(int i=0;i<obj_size;i++)
    {
        float t1=INFINITY, t2 =INFINITY;
        if(objects[i]->intersect(ray, t1, t2)){
            if(t1<0) t1=t2;
            if(t1<t_min){
                t_min = t1;
                object = objects[i];
            }
        }
    }
    for(int i=0;i<light_size;i++){
        float t1=INFINITY, t2=INFINITY;
        if(lightings[i]->intersect(ray, t1, t2)){
            if(t1<0) t1=t2;
            if(t1<t_min){
                t_min = t1;
                object = lightings[i];
            }
        }
    }
    
    float t = INFINITY;
    int face_id =-1;
    if(model) face_id = model->intersect(ray.origin, ray.direction, t);
    if (t < t_min) {
        t_min = t;
        object = NULL;
    }
    
    Vec3f surfColor = 0;
    
    if (face_id >= 0) {
        surfColor = model->surfColor * model->faceNorms[face_id].dot(-ray.direction);
        return surfColor;
    }
    
    //光线与场景中的物体有相交
    if(object != NULL){
        
        Vec3f phit = ray.getPoint(t_min);
        Vec3f nhit = object->nhit(phit);
        
        float IdotN = ray.direction.dot(nhit);
        float facingratio = std::max(float(0), -IdotN);
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.3);//菲涅尔效应
        
        float rayProjection = abs(ray.direction.dot(nhit));
        Vec3f reflDirection = nhit * 2*rayProjection + ray.direction;
        Ray reflectRay = Ray(phit,reflDirection);
        
        color +=  object->getTextureColor(phit) *object->surfaceColor  *object->Ka * AmbientLightColor *AmbientStrength;
        
        color += object->emissionColor;
        
        if(object->Ksp>0 && depth<MAX_DEPTH){
            
            //发生镜面反射，Phong光照模型，光线追踪
            Vec3f reflColor = trace(reflectRay,objects,lightings,model,depth+1);
            
            Vec3f h = (reflDirection - ray.direction).normalize();
            
            color += ( reflColor* fresneleffect  *std::max(0.0, pow(h.dot(nhit),object->Fsp))) *object->Ksp ;//镜面反射和surfacecolor无关
        }
        
        if(object->Kd>0){
            double shadow =1.0;
            for(int i=0;i<light_size;i++){
                float D = (lightings[i]->center-phit).length()-lightings[i]->radius;
                float atten = 1/(1+0.03*D + 0.001 * D *D);//光强随距离递减
                Vec3f s = (lightings[i]->center-phit).normalize();
                Vec3f origin = phit + (nhit * bias);
                Ray lightRay(origin,s);
                Vec3f transmission(1.0f);
                Object *shadowCast = NULL;
                float tnear = INFINITY;
                
                for(int j =0;j<obj_size;j++){
                    float t1=INFINITY, t2=INFINITY;
                    if(objects[j]->intersect(lightRay, t1, t2)){
                        if(t1<tnear){
                            tnear = t1;
                            shadowCast = objects[j];
                        }
                    }
                }
                for(int j =0;j<light_size;j++){
                    float t1=INFINITY, t2=INFINITY;
                    if(i!=j && lightings[j]->intersect(lightRay, t1, t2)){
                        if(t1<tnear){
                            tnear = t1;
                            shadowCast = lightings[j];
                        }
                    }
                }
                
                if(shadowCast){
                    shadow = std::max(0.0, shadow - (1.0 - shadowCast->transparency));
                    transmission = transmission * shadowCast->surfaceColor * shadow;
                }
                if(model){
                float tmp_t;
                int tmp_face_id = model->intersect(lightRay.origin,lightRay.direction, tmp_t);
                if (tmp_face_id > 0) transmission = 0;
                }
                
                color +=  object->getTextureColor(phit) *  object->surfaceColor * lightings[i]->emissionColor *transmission* atten * std::max(0.0f, s.dot(nhit)) *object->Kd;
            }
        }
        
    }
    
    return color;
}

void render(const std::vector<Object *> objects, const std::vector<Sphere *>lightings, Model *model,Camera &camera, GLFWwindow *window){
    
    
    unsigned char* pix = new unsigned char[WINDOW_WIDTH*WINDOW_HEIGHT * 3];
    
    Vec3f *image = new Vec3f[WINDOW_WIDTH * WINDOW_HEIGHT], *pixel = image;
    
    for (unsigned y = 0; y < WINDOW_HEIGHT; ++y) {
        for (unsigned x = 0; x < WINDOW_WIDTH; ++x, ++pixel) {
            Ray ray=camera.generateRay(x, y, WINDOW_WIDTH, WINDOW_HEIGHT);
            int depth =0;
            Vec3f color = trace(ray, objects,lightings,model,depth);
            *pixel = color;
        }
    }
    
    pixel = image;
    for (unsigned i = 0; i < WINDOW_HEIGHT; i++)
        for (unsigned j = 0; j < WINDOW_WIDTH; j++, pixel++){
            pix[3* (i*WINDOW_WIDTH + j)] = std::min(pixel->x, float(1)) * 255;
            pix[3 * (i*WINDOW_WIDTH + j) + 1] = std::min(pixel->y, float(1)) * 255;
            pix[3 * (i*WINDOW_WIDTH + j) + 2] = std::min(pixel->z, float(1)) * 255;
        }
    
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawPixels(WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pix);
        glfwSwapBuffers(window);
    }
    delete[] pix;
    glfwTerminate();
}

int main(int argc, const char * argv[]) {
    
    std::vector<Object *> objects;
    std::vector<Sphere *>lightings;
    
    char imagePath[200];
    char modelPath[200];
    char defaultimagePath[200] = "/Users/sunchenge/Desktop/计算机图形学/raytracing/raytracing/1.bmp";
    char defaultmodelPath[200] = "/Users/sunchenge/Desktop/计算机图形学/raytracing/raytracing/bunny_100.txt";
    
    std::cout<<"input texture image path:";
    std::cin>>imagePath;
    std::cout<<"input model path:";
    std::cin>>modelPath;
    
    if(strlen(imagePath)<2){//没什么用，自己运行时偷懒用的
        strcpy(imagePath,defaultimagePath);
        strcpy(modelPath,defaultmodelPath);
    }
    
    
    
    lightings.push_back(new Sphere(Vec3f(20,12,20),5,Vec3f(1.0f,1.0f,1.0f),1,0,WHITE,0.3,1.0,0.0));

    
    Cube * cube1 = new Cube(Vec3f(4,9,0),3,3,3,Vec3f(1.0f,1.0f,1.0f),1,0);
    Cube * cube2 = new Cube(Vec3f(7,0.5,0),1.5,1.5,1.5,Vec3f(1.0f,0.0f,0.0f),1,0);
    cube1->setTextureImage(imagePath);
    objects.push_back(cube1);
    objects.push_back(cube2);
    Sphere * sphere1 =new Sphere(Vec3f(-7,8,3),3,Vec3f(1.0f,1.0f,1.0f),1,0);
    sphere1->setTextureMapping(CHESS);
    objects.push_back(sphere1);
    objects.push_back(new Sphere(Vec3f(7,4,1.5),1.5,Vec3f(0.5f,1.0f,1.0f),1,0));
    objects.push_back(new Cube(Vec3f(4.5,9.3,3),2,2,3,Vec3f(1.0f,0.0f,0.0f),1,0));
    objects.push_back(new Plane(Vec3f(0.0f,0.0f,1.0f),0.0f,Vec3f(0.8f,0.8f,0.8f),1,0));
    Model *bunny = new Model(modelPath);

    
    Camera camera(Vec3f(20,5,5), Vec3f(5,5,5), Vec3f(-1.73,0,1),60);
    
    glewExperimental=true;
    if( !glfwInit())
    {
        return -1;
    }
//    glfwWindowHint(GLFW_SAMPLES, 4);
//    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
//    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
//    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
//    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH/2, WINDOW_HEIGHT/2, "Ray Tracing", nullptr, nullptr);//玄学 不除2就只会在左下角显示
    if (window == nullptr)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    

    glfwMakeContextCurrent(window);
    
    glewExperimental=true;
    if(glewInit()!=GLEW_OK)
    {
        return -1;
    }
    
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    
    render(objects,lightings,bunny,camera,window);
    
    return 0;
    
}
