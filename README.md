# Ray Tracing
* Assignment2 of CS337 in SJTU
* By Chenge Sun
* With C++ glew+glfw3 Xcode(ver10.1) on macOS Mojave(ver10.14.2)

## Including Files
* In executablefile :
	* *raytracing* (an UNIX executable file)
* Project file *raytracing.xcodeproj*
* In raytracing file:
	* *bmpstruct.h* used to read *.bmp* files
	* *GLVector.h* creates a 3D vector data structure for convenience
	* *GLObject.h* creates :
		* Camera
		* Light Ray
		* Objects: Sphere,Plane,Cube,Model
	* *main.cpp* implements tracing and rendering, creates a scene and displays it.
	* *bunny\_100.txt*: the vertice and faces information of Stanford Bunny Model with 100 faces
	* *1.bmp* a demo image texture
* In image file:
	* screenshots of the demo scene

## How to Run
### Xcode
1. install glew+glfw3
2. project -> Target -> Build Setting 
	
	Add `/usr/local/include` to `header search paths`
3. Build Phases

	Add `opengl.framework,libGLEW.dylib,libGlfw.dylib` to `Link binary with libraries`
4. input a texture image path and a model path for my demo scene. if you want to use my demo image and model, input using the following form:

~~~
/path/to/this/file/raytracing/1.bmp
/path/to/this/file/raytracing/bunny_100.txt
~~~

Then you can see the demo scene as shown in */images/4\_Scene* after a time-consuming rendering process

### Terminal
Open Terminal on MacOS and do the following instructions:

~~~
cd /path/to/this/file/executablefile
./raytracing
input texture image path:/path/to/this/file/raytracing/1.bmp
input model path:/path/to/this/file/raytracing/bunny_100.txt
~~~

Then you can see the demo scene as shown in */images/4\_Scene* after a time-consuming rendering process
		
## Finished Work
* Ray trace and draw a plane, some sphere objects, some cube objects and a Stanford Bunny model object.
* Display a point lighting using a glowing sphere
* Perform Phong-shading, display shadows and recursive reflection(up to 3 interations)
* Read in material properties of sphere and cube objects including:
	* Surface color
	* Emission color
	* Reflection coefficient( $k_a,k_d,k_{sp},sunshiness$ )
	* Image files for implementing texture mapping(only support *.bmp* files currently)
	* Build-in texture（only support chessboard texture currently）
* Read in model files(in *.txt* files with vertices and faces of the model)

## Acknowledgement
I learnt how to read in *.bmp* file on MacOS from [here](https://github.com/hggq/c-bmp-class).

The model of Stanford Bunny is from [here](http://www.cc.gatech.edu/projects/large_models/bunny.html)

## Contact
sunchenge1997@sjtu.edu.cn
