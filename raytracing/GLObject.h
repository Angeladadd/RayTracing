//
//  GLObject.h
//  raytracing
//
//  Created by sunchenge on 2018/12/28.
//  Copyright © 2018 sunchenge. All rights reserved.
//

#ifndef GLObject_h
#define GLObject_h

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
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include "bmpstruct.h"

#define PI 3.141592653589
#define   WIDTHBYTES(bits) (((bits)+31)/32*4)
#define smallEnough 1e-4



using namespace bmp;

Vec3f WHITE(1);
Vec3f BLACK(0);

//以后可添加
enum TEXTURE{
    CHESS
};


//光线类
class Ray{
public:
    Vec3f origin;//光源的位置
    Vec3f direction;//光线的方向向量，已被标准化
    float intense;
    Ray(Vec3f &o, Vec3f &d, float i=1.0f):origin(o),intense(i){
        direction = d.normalize();
    }
    Vec3f getPoint(float t){
        return origin + direction * t;
    }
};

    //相机
class Camera{
public:
    Vec3f eye;
    Vec3f center;
    Vec3f up;
    float FOV;
    Camera(Vec3f e, Vec3f c, Vec3f u, int fov):eye(e),center(c),up(u), FOV(fov/180.0*PI){
    }
    //PPT光线算法
    Ray generateRay(int i, int j, int width, int height){
        Vec3f a = eye - center;
        Vec3f w = a.normalize();
        Vec3f b = up;
        Vec3f u = b.cross(w);
        u = u.normalize();
        Vec3f v = w.cross(u);
        float alpha = tan(FOV/2)*(i-width/2)/(width/2);
        float beta = tan(FOV/2)*(j-height/2)/(height/2);
        Vec3f direction = u*alpha+v*beta-w;
        Ray ray(eye, direction);
        return ray;
    }
};



class Object {
public:
    Vec3f surfaceColor;
    Vec3f emissionColor;//光源用
    float reflection;//反射系数，后来被抛弃了，但没时间删了
    float transparency;//目前没来得及实现透明物体
    float Ka;//环境光反射系数
    float Kd;//漫反射系数
    float Ksp;//镜面反射系数
    int Fsp;//sunshineness镜面反射指数
    bool textureMap;//是否有纹理映射
    enum TEXTURE texture;//内置纹理映射的种类
    Object(Vec3f sc,float refl, float trans,Vec3f ec=BLACK,float ka=0.3, float kd=0.8, float ksp=0.5, int fsp=6):surfaceColor(sc),emissionColor(ec),reflection(refl),transparency(trans), Ka(ka),Kd(kd),Ksp(ksp),Fsp(fsp),textureMap(false){}
    virtual bool intersect(Ray &ray, float &t1, float &t2) const = 0;
    virtual Vec3f nhit(Vec3f &phit) const = 0;
    virtual void setTextureImage(char * imageName) = 0;
    virtual void setTextureMapping(enum TEXTURE tex) = 0;
    virtual Vec3f getTextureColor(Vec3f &phit) const=0;
};


class Sphere:public Object{
public:
    Vec3f center;
    float radius;
    float radius2;
    Sphere(Vec3f c, float r, Vec3f sc, float refl, float trans, Vec3f ec=BLACK,float ka=0.3, float kd=1.0, float ksp=0.9, int fsp=50):center(c),radius(r),radius2(r*r),Object(sc,refl,trans,ec,ka,kd,ksp,fsp){}
    bool intersect(Ray &ray, float &t1, float &t2) const{
        Vec3f distance = center - ray.origin;
        float distanceProjection = distance.dot(ray.direction);
        if(distanceProjection<0) return false;//不相交
        float distanceLength2 = distance.length2();
        float rdistance2 = distanceLength2-distanceProjection*distanceProjection;
        if(rdistance2>radius2) return false;//不相交
        float interception = sqrt(radius2 - rdistance2);//截距的一半
        t1 = distanceProjection - interception;
        t2 = distanceProjection + interception;
        return true;
    }
    //标准化过的法线
    Vec3f nhit(Vec3f &phit) const {
        return (phit - center).normalize();
    }
    void setTextureImage(char * imageName) {
        textureMap = 1;
        //不想写了
    }
    void setTextureMapping(enum TEXTURE tex) {
        textureMap = 1;
        texture = tex;
    }
    Vec3f getTextureColor(Vec3f &phit) const{
        if(!textureMap) return WHITE;
        Vec3f sphereCord;
        Vec3f relativeCord = phit-center;
        sphereCord.x = radius;
        sphereCord.y = PI-acos(relativeCord.z/radius);//0-PI
        sphereCord.z = atan(relativeCord.y/relativeCord.x);//-PI/2-PI/2,下面范围讨论后变成0-2PI
        float yy=relativeCord.dot(yVec);
        float xx=relativeCord.dot(xVec);
        if(yy<0 && xx<0) sphereCord.z += PI;
        else {
            if(yy<0) sphereCord.z += 2*PI;
            else if(xx<0) sphereCord.z += PI;
        }
        
        float scale = PI/8;
        
        bool flag = 0;
        for(int i=0;i<7;i+=2){
            if(sphereCord.y<i*scale) break;
            if(sphereCord.y>i*scale && sphereCord.y<(i+1)*scale){
                flag = 1;
                break;
            }
        }
        for(int i=0;i<15;i+=2){
            if(sphereCord.z<i*scale) break;
            if(sphereCord.z>i*scale && sphereCord.z<(i+1)*scale){
                flag = !flag;//如果两个同时黑色就变成白色 模拟棋盘
                break;
            }
        }
        if(flag) return BLACK;
        return WHITE;
    }
};

class Plane:public Object{
public:
    Vec3f normal;//平面法向量 判断朝向
    float distance;//与原点距离
    Plane(Vec3f n, float d, Vec3f sc,  float refl, float trans, Vec3f ec=BLACK,float ka=0.3, float kd=1.0, float ksp=0.5f, int fsp=6):normal(n.normalize()), distance(d),Object(sc,refl,trans,ec,ka,kd,ksp,fsp){}
    bool intersect(Ray &ray, float &t1, float &t2) const{
        if(normal.dot(ray.direction)>=0) return false;//相交点积小于0
        float o_d = abs(ray.origin.dot(normal) - distance);
        float cosa = - ray.direction.dot(normal);
        t1 = o_d/cosa;
        t2 = t1;
        return true;
    }
    void setTextureImage(char * imageName) {
        //不想写了
    }
    void setTextureMapping(enum TEXTURE tex){
        //不想写了
    }
    Vec3f getTextureColor(Vec3f &phit) const{
        return WHITE;
    }
    Vec3f nhit(Vec3f &phit) const{
        return normal;
    }
};

class Cube:public Object{
public:
    Vec3f vertex;
    float xLength;
    float yLength;
    float zLength;
    
    int imageW;//纹理图片的信息
    int imageH;
    Vec3f *imageBuffer;//颜色
    
    Cube(Vec3f v, float xl, float yl, float zl, Vec3f sc,  float refl, float trans, Vec3f ec=BLACK,float ka=0.3, float kd=1.0, float ksp=0.9, int fsp=6):vertex(v),xLength(xl),yLength(yl),zLength(zl),Object(sc,refl,trans,ec,ka,kd,ksp,fsp){}
    void setTextureImage(char * imageName) {
            //还要读BMP文件 这个作业好麻烦！！！！！
        textureMap = 1;
        bmp::bitbmp mybmp(imageName);
        mybmp.read(imageName);
        imageW = mybmp.width;
        imageH = mybmp.height;
        
        imageBuffer = new Vec3f[imageW * imageH];
        for(int i=0;i<imageH;i++)
            for(int j=0;j<imageW;j++){
                imageBuffer[i*imageW + j].x = (mybmp.getData()[i*imageW +j].Red-'0')/255.0f;//bmp中是0-255 把它标准化
                imageBuffer[i*imageW + j].y = (mybmp.getData()[i*imageW +j].Green-'0')/255.0f;
                imageBuffer[i*imageW + j].z = (mybmp.getData()[i*imageW +j].Blue-'0')/255.0f;
            }

    }
    void setTextureMapping(enum TEXTURE tex) {
        textureMap = 1;
        
    }
    Vec3f getTextureColor(Vec3f &phit) const{
        if(!textureMap) return WHITE;
        
        float x=0,y=0;
        float top ;
        float bottom ;
        float left ;
        float right ;
        
        if((abs(phit.x-vertex.x-xLength)<smallEnough || abs(phit.x-vertex.x)<smallEnough)){
            x = phit.y;
            y = phit.z;
            top = vertex.z+zLength;
            bottom = vertex.z;
            left = vertex.y;
            right = vertex.y+yLength;
        }
        else {
            if((abs(phit.y-vertex.y-yLength)<smallEnough || abs(phit.y-vertex.y)<smallEnough)){
                x = phit.x;
                y = phit.z;
                top = vertex.z+zLength;
                bottom = vertex.z;
                left = vertex.x;
                right = vertex.x+xLength;
            }
            else {
                if((abs(phit.z-vertex.z-zLength)<smallEnough || abs(phit.z-vertex.z)<smallEnough)){
                    x = phit.x;
                    y = phit.y;
                    top = vertex.y+yLength;
                    bottom = vertex.y;
                    left = vertex.x;
                    right = vertex.x+xLength;
                }
            else return WHITE;
            }
        }
            
        
        
        float u = (x-left)/(right-left);
        float v = (y-bottom)/(top-bottom);
        
        int intu = (int)imageW * u;//索引里面是整数所以这里有近似，如果用颜色很精细的图会变形，大面积同样颜色没什么问题
        int intv = (int)imageH * v;
        
        return imageBuffer[intv*imageW + intu];
        
    }

    bool intersect(Ray &ray, float &t1, float &t2) const {
        //判断外部相交
        
        Vec3f nTop(0.0f,0.0f,1.0f);//正方体六个面法向量
        Vec3f nBottom(0.0f,0.0f,-1.0f);
        Vec3f nLeft(0.0f,-1.0f,0.0f);
        Vec3f nRight(0.0f,1.0f,0.0f);
        Vec3f nFront(1.0f,0.0f,0.0f);
        Vec3f nBack(-1.0f,0.0f,0.0f);
        
        std::vector<Vec3f> normals;
        normals.push_back(nTop);
        normals.push_back(nBottom);
        normals.push_back(nLeft);
        normals.push_back(nRight);
        normals.push_back(nFront);
        normals.push_back(nBack);
        
        for(int i=0;i<6;i++){
            if(ray.direction.dot(normals[i])>= 0) continue;
            float cosa = abs(ray.direction.dot(normals[i]));
            float t;
            bool flag = 0;
            Vec3f tmp_p;
            switch(i){
                case 0:
                    t = (ray.origin.z-vertex.z-zLength)/cosa;//这里还要判断光线的origin在平面的哪个方向
                    if(t<0) break;
                    tmp_p = ray.getPoint(t);
                    flag = tmp_p.x>=vertex.x && tmp_p.x<=vertex.x+xLength && tmp_p.y>=vertex.y && tmp_p.y<=vertex.y+yLength;//另外两个坐标的范围是不是在正方形里
                    break;
                case 1:
                    t = (vertex.z-ray.origin.z)/cosa;
                    if(t<0) break;
                    tmp_p = ray.getPoint(t);
                    flag = tmp_p.x>=vertex.x && tmp_p.x<=vertex.x+xLength && tmp_p.y>=vertex.y && tmp_p.y<=vertex.y+yLength;
                    break;
                case 2:
                    t = (vertex.y-ray.origin.y)/cosa;
                    if(t<0) break;
                    tmp_p = ray.getPoint(t);
                    flag = tmp_p.z>=vertex.z && tmp_p.z<=vertex.z+zLength && tmp_p.x>=vertex.x && tmp_p.x<=vertex.x+xLength;
                    break;
                case 3:
                    t = (ray.origin.y-vertex.y-yLength)/cosa;
                    if(t<0) break;
                    tmp_p = ray.getPoint(t);
                    flag = tmp_p.z>=vertex.z && tmp_p.z<=vertex.z+zLength && tmp_p.x>=vertex.x && tmp_p.x<=vertex.x+xLength;
                    break;
                case 4:
                    t = (ray.origin.x-vertex.x-xLength)/cosa;
                    if(t<0) break;
                    tmp_p = ray.getPoint(t);
                    flag = tmp_p.z>=vertex.z && tmp_p.z<=vertex.z+zLength && tmp_p.y>=vertex.y && tmp_p.y<=vertex.y+yLength;
                    break;
                default:
                    t = (vertex.x-ray.origin.x)/cosa;
                    if(t<0) break;
                    tmp_p = ray.getPoint(t);
                    flag = tmp_p.z>=vertex.z && tmp_p.z<=vertex.z+zLength && tmp_p.y>=vertex.y && tmp_p.y<=vertex.y+yLength;
                    break;
            }
            if(!flag) continue;
            if(t<t1){
                t2 = t1;
                t1=t;
            }
            else t2 = t;
        }
        
        if(t1<INFINITY) return true;
        else return false;
    }
    
    Vec3f nhit(Vec3f &phit) const{
        if(abs(phit.x-vertex.x)<smallEnough) return Vec3f(-1.0f,0.0f,0.0f);//不精确直接等会出问题
        if(abs(phit.y - vertex.y)<smallEnough) return Vec3f(0.0f,-1.0f,0.0f);
        if(abs(phit.z - vertex.z)<smallEnough) return Vec3f(0.0f,0.0f,-1.0f);
        if(abs(phit.x-vertex.x - xLength)<smallEnough) return Vec3f(1.0f,0.0f,0.0f);
        if(abs(phit.y-vertex.y - yLength)<smallEnough) return Vec3f(0.0f,1.0f,0.0f);
        else return Vec3f(0.0f,0.0f,1.0f);

    }
        
};
class Model {
public:
    int nVer, nFace;                //顶点数和片元数
    vector<Vec3f> vers;                //顶点
    vector<Vec3d> faces;            //片元对应的顶点索引
    vector<Vec3f> faceNorms;        //片元法向量
    Vec3f surfColor;
    
    Model(char* filename) {
        ifstream fin(filename);
        if(!fin) {std::cout<<"open failed!\n"; exit(0);}
        char tmpchar[100];
        int tmpnum;
        surfColor.x = 1;  surfColor.y = 1; surfColor.z = 0.7;
        fin >> tmpchar;
        fin >> nVer >> nFace >> tmpnum;
        for (int i = 0; i < nVer; ++i) {
            Vec3f pt;
            fin >> pt.y >> pt.z >> pt.x;//文件的坐标方向和我是反的
            pt = pt*30 + Vec3f(-3,-3,-1);
            vers.push_back(pt);
        }
        for (int i = 0; i < nFace; ++i) {
            Vec3d fs;
            Vec3f nm;
            fin >> tmpnum >> fs.x >> fs.y >> fs.z;
            faces.push_back(fs);
            Vec3f v1 = vers[fs.y] - vers[fs.x], v2 = vers[fs.z] - vers[fs.y];
            nm = v1.cross(v2);
            nm.normalize();
            faceNorms.push_back(nm);
        }
        fin.close();
    }
    
    int intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t) {
        for (int i = 0; i < nFace; ++i) {
            if (intersect(i, rayorig, raydir, t))
                return i;
        }
        return -1;
    }
    
private:
    bool intersect(int i, const Vec3f &rs, const Vec3f &rd, float &t) {
        Vec3f nm = faceNorms[i], pt1 = vers[faces[i].x], pt2 = vers[faces[i].y], pt3 = vers[faces[i].z];
        float tmp_t;
        float dir_dot_norm = rd.dot(nm);
        if (dir_dot_norm >= 0) {
            return false;
        }
        tmp_t = nm.dot(pt1 - rs) * (1 / dir_dot_norm);
        
        Vec3f P = rs + rd * tmp_t;
        //判断是否在三角形里
        
        float areaAll = ((pt2 - pt1).cross(pt3 - pt2)).length();
        float area1 = ((pt2 - pt1).cross(P - pt2)).length();
        float area2 = ((pt3 - pt2).cross(P - pt3)).length();
        float area3 = ((pt1 - pt3).cross(P - pt1)).length();
        if (area1 + area2 + area3 > areaAll + area1 / 1000) {

            return false;
        }
        
        t = tmp_t;
        return true;
    }
};


#endif /* GLObject_h */

