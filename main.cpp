#include <iostream>
#include "src\core\vec.hpp"
#include "src\core\camera.hpp"
#include "src\core\scene.hpp"
#include "src\core\integrator.hpp"
#include "src\core\light.hpp"

#include "src\primitives\sphere.hpp"
#include "src\primitives\box.hpp"
#include "src\primitives\plane.hpp"
#include "src\primitives\disk.hpp"
#include "src\primitives\triangle.hpp"
#include "src\primitives\parallelogram.hpp"
#include "src\primitives\csg.hpp"
#include "src\primitives\quadric.hpp"



#include "external/tigr/tigr.h"
#ifdef _MSC_VER
#pragma comment(lib,"D3D9.lib")
#endif
using namespace HydrogenCG;

/*

	DONE:
	- add animations and physics
	- quadrics morphing
	- quadrics
	- collisions
	- non-uniform scale
	- texture image glsl like
	- vec2
	- bsdf
	- fixed CSG
	- added distance function
	- fixed instance
	- fixed intersection structure
	- Do not nest instances, permitted only for composite objects - CSG, plain objects needs to remain plain
	- otherwise you nest bsdfs and transforms
	- collection pivot?

	NOT IMPLEMENTED:
	- fog integrator (scrapped, possibly in the raymarching section_
	- fake subsurface scattering integrator (scrapped, possibly in the ray marching section)
	- area light (postponed till distributed rt section)
	- AO (postponed till pt section)
	- bicubic filtering and other filters (left as an exercise for the reader)
	- cubics, quartics (torus) (scrapped, left as an exercise for the reader - add notes on tensor contractions, link to classification)
	- other surfaces (patches, subdivision) (scrapped entirely)
	- 3d mesh loading and rendering (idk?)
	- acceleration structures (possibly later, linked to 3d mesh loading and rendering really)
	- bent rays (quadratic, cubic) (seems fun, but idk)
	- different camera models (postponed till distributed rt section, possibly here?)
	- other brdfs (will introduce some in the rt and some in the pt section)
	- environment maps (idk?)
	- samplers (distributed rt or pt section)
	- signal reconstrution (distributed rt or pt section)
	- adaptive sampling (pt section?)
*/

#include <chrono>

class Stopwatch
{
private:
	std::chrono::high_resolution_clock clock;
	std::chrono::high_resolution_clock::time_point timeStamp;

public:
	void reset() { timeStamp = clock.now(); }
	float elapsed() 
	{
		std::chrono::duration<float> dur = clock.now() - timeStamp;
		return dur.count();
	}
};

int main()
{
	// set up timers
	Stopwatch fpsWatch;		// timer used for fps measurement
	Stopwatch deltaWatch;	// timer used for frame time measurement
	Stopwatch time;			// timer used for total elapsed time measurement

	// set up the image that we will render into
	Image image;
	image.init(320,240);

	// set up 
	Camera camera;
	Scene scene;


	// create a checkerboard texture
	auto checker = bsdfDiffuse(textureCheckerboard(vec3(1), vec3(0)));

	// create a "solar system" transformation hierarchy

	// add a pivot for the earth
	auto earth = collection();
	// add a sphere representing the Earth
	earth->add(sphere(vec3(0), 0.1f, vec3(0,0.3,1.0f)));
	// add a sphere representing the moon
	earth->add(sphere(vec3(0.2f,0,0), 0.05f, vec3(0.6f)));
	// translate the pivot for the Earth (both the earth and the moon move through this action)
	earth->translate(vec3(0, 0, 3));

	// create a Sun pivot
	auto sun = collection();
	sun->translate(vec3(-3, 3, 5));
	// create a Sphere representing the sun
	sun->add(sphere(vec3(0), 1.0f, vec3(1,1,0)));
	// attach the Earth hierarchy to the Sun pivot
	sun->add(earth);
	
	// Create a box and a sphere
	auto tempBox = box({ -1,-1,-1 }, { 1,1,1 }, { 1, 0.3f, 0 });
	auto tempSphere = sphere({ 0,0,0 }, 1.3f, { 0, 0, 1 });
	// carve out the sphere from the box, then translate the result
	auto carvedOutBox = tempBox - tempSphere + vec3(0, 0, 6);
	// create a plane for the ground, and apply the checker texture to it
	auto groundPlane = plane(vec3(0, -2, 0),vec3(0,0,1), normalize(vec3(1,0,1)), checker);

	// create the morphing quadric
	/*auto morphingQuadric = quadricHyperbolicCylinder();
	morphingQuadric->bsdf = checker;*/
	
	mat4 initialMorphingQuadricMatrix = mat4::identity;
	initialMorphingQuadricMatrix.at(3, 3) = -1;
	auto morphingQuadricSurface = std::make_shared<Quadric>(vec3(0,0,0), initialMorphingQuadricMatrix);

	// create a sequence of quadric matrices into which the quadric will morph into (cyclically)
	std::vector<mat4> quadricMorphs;
	int currentMorph = 0;
	// A sphere 
	quadricMorphs.push_back(mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1)));
	// A cylinder 
	quadricMorphs.push_back(mat4(vec4(1, 0, 0, 0), vec4(0, 0, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1)));
	// A cylinder 
	quadricMorphs.push_back(mat4(vec4(1, 0, 0, 0), vec4(0, 0, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1)));
	// A one sheet hyperboloid
	quadricMorphs.push_back(mat4(vec4(1, 0, 0, 0), vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1)));
	// A cone
	quadricMorphs.push_back(mat4(vec4(1, 0, 0, 0), vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, 0)));
	// A two sheet hyperboloid
	quadricMorphs.push_back(mat4(vec4(1, 0, 0, 0), vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, 1)));
	// Two crossed planes
	quadricMorphs.push_back(mat4(vec4(0, 0, 0, 0), vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, 0)));
	// A hyperbolic cylinder
	quadricMorphs.push_back(mat4(vec4(0, 0, 0, 0), vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0,0,0, -1)));
	// A hyperbolic paraboloid
	quadricMorphs.push_back(mat4(vec4(0, 0, 0, -1), vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(-1, 0, 0, 0)));
	// A parabolic cylinder
	quadricMorphs.push_back(mat4(vec4(0, 0, 0, -1), vec4(0, 0, 0, 0), vec4(0, 0, 1, 0), vec4(-1, 0, 0, 0)));
	// A paraboloid
	quadricMorphs.push_back(mat4(vec4(0, 0, 0, -1), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(-1, 0, 0, 0)));

	// get a reference to the matrix of the morphing quadric
	mat4& morphingQuadricMatrix = morphingQuadricSurface->m;

	// add the ground and carved out box to the scene
	scene.collection.add(groundPlane);
	scene.collection.add(carvedOutBox);
	scene.collection.add(sun);
	
	auto morphingQuadric = std::make_shared<SurfaceInstanceBSDF>(morphingQuadricSurface, checker);
	// create the morphing quadric instance contained withing a box
	scene.collection.add((1.0f*morphingQuadric) * (3.0f*tempBox) + vec3(5,3,2));
	
	// add a point light to the scene
	scene.lights.push_back(std::make_unique<LightPoint>(vec3( 0, 3, 0 ),100* vec3(1.0f)));

	// put all integrators in a vector
	std::vector<std::unique_ptr<Integrator>> integrators;
	integrators.push_back(std::make_unique<Integrator>());
	integrators.push_back(std::make_unique<IntegratorDepth>());
	integrators.push_back(std::make_unique<IntegratorNormal>());
	integrators.push_back(std::make_unique<IntegratorColor>());
	integrators.push_back(std::make_unique<IntegratorTransparency>());
	integrators.push_back(std::make_unique<IntegratorDiffuseLocal>());
	integrators.push_back(std::make_unique<IntegratorDiffuse>());
	int currentIntegratorIndex = 0;

	Tigr *screen = tigrWindow(image.w(), image.h(), "HydrogenCG", 0);
	Tigr* tigrBmp = tigrBitmap(image.w(), image.h());
	
	

	// save old mouse position and button state for drag rotation
	double oldMouseX = 0;
	double oldMouseY = 0;
	bool mouseButtonDownPrevFrame = false;
	// save azimuthal and polar angle for camera rotation
	float azimuthalAngle = 0;
	float polarAngle = 0;
	float oldAzimuthalAngle = 0;
	float oldPolarAngle = 0;
	
	// time and fps
	u32 frames = 0;
	u32 fps = 0;
	fpsWatch.reset();
	time.reset();
	// set the first morph for the morphing quadric
	float morphLastT = time.elapsed();
	while (!tigrClosed(screen) && !tigrKeyDown(screen, TK_ESCAPE))
	{
		// save the time at the beginning of the frame to measure delta time
		deltaWatch.reset();

		// measure fps
		if (fpsWatch.elapsed() > 1.0f)
		{
			fps = frames;
			fpsWatch.reset();
			frames = 0;
		}

		// animate
		float currT = sinf(time.elapsed());
		scene.collection[1]->rotate(vec3(1.0f*currT, 2.0f*currT, 0));
		scene.collection[1]->translate(vec3(0, 0.1f*currT, 0));
		scene.collection[1]->scale(0.3f+smoothstepQuintic(fabs(currT)));
		sun->rotate(vec3(0, time.elapsed(), 0));
		earth->rotate(vec3(0, 10.0f*time.elapsed(), 0));
		
		// morph quadric
		float morphNewT = 1.0f*time.elapsed();
		float t = morphNewT - morphLastT;
		if (t > 1.0f)
		{
			morphLastT = morphNewT;
			currentMorph = (currentMorph + 1) % quadricMorphs.size();
		}
		else
		{
			int nextMorph = (currentMorph + 1) % quadricMorphs.size();
			morphingQuadricMatrix = (1 - t)*quadricMorphs[currentMorph] + t * quadricMorphs[nextMorph];
		}
		
		
		// render with the current integrator
		integrators[currentIntegratorIndex]->render(image, camera, scene);

		// gamma correct and write image to buffer
		for (u32 i = 0; i < image.w()*image.h(); ++i)
		{
			tigrBmp->pix[i].a = 255;
			vec3 col = image[i];
			// throw away negative values
			col = max(col, 0.0f);
			// apply gamma
			col = pow(col, 1.0f/2.2f);
			// map from [0,1] to [0,maxVal] and round
			col = round(float(255)*col);
			// throw away values higher than maxVal
			col = min(col, float(255));
			u32 r = static_cast<u32>(col.r);
			u32 g = static_cast<u32>(col.g);
			u32 b = static_cast<u32>(col.b);
			tigrBmp->pix[i].r = r;
			tigrBmp->pix[i].g = g;
			tigrBmp->pix[i].b = b;
		}
		tigrBlit(screen, tigrBmp, 0, 0, 0, 0, image.w(), image.h());

		// output fps
		tigrPrint(screen, tfont, 0, 0, tigrRGB(255, 255, 255), "Fps: %d", fps);
		tigrUpdate(screen);

		// measure delta time
		float delta = deltaWatch.elapsed();

		// camera speed
		vec3 speed = 5.0f*vec3(1, 1, 1);

		// compute the closest distance from the camera to the scene
		float sceneDist = scene.distance(camera.o);
		// estimate (numerically, using finite differences) the normal from the camera to the scene
		vec3 norm = normalize(vec3(scene.distance(camera.o + vec3(EPSILON, 0, 0)) - scene.distance(camera.o - vec3(EPSILON, 0, 0)),
			scene.distance(camera.o + vec3(0, EPSILON, 0)) - scene.distance(camera.o - vec3(0, EPSILON, 0)),
			scene.distance(camera.o + vec3(0,0,EPSILON)) - scene.distance(camera.o - vec3(0,0,EPSILON))));
		// if the camera has entered the scene - push it back along the estimated normal
		camera.o += max(0.2f-sceneDist,0.0f)*norm;

		// camera movement
		if (tigrKeyHeld(screen, 'W')) camera.o += delta*speed.z*camera.e2; // forward
		if (tigrKeyHeld(screen, 'S')) camera.o += -delta*speed.z*camera.e2; // backward
		if (tigrKeyHeld(screen, 'A')) camera.o += -delta*speed.x*camera.e0; // strafe left
		if (tigrKeyHeld(screen, 'D')) camera.o += delta*speed.x*camera.e0;  // strafe right
		if(tigrKeyHeld(screen, 'Q')) camera.o += delta * speed.y*camera.e1; // up
		if (tigrKeyHeld(screen, 'E')) camera.o += -delta * speed.y*camera.e1; // down
		// change rendering mode
		if (tigrKeyDown(screen, TK_SPACE)) currentIntegratorIndex = (currentIntegratorIndex + 1) % integrators.size();

		// mouse input
		int inputMouseX, inputMouseY, inputMouseButtons;
		tigrMouse(screen, &inputMouseX, &inputMouseY, &inputMouseButtons);
		// map mouse coordinates from [0,width] x [0,height] to [0,1]^2
		double newMouseX = 2.0*inputMouseX / double(image.w()) - 1.0;
		double newMouseY = -2.0*inputMouseY / double(image.h()) + 1.0;
		float diffMouseX = 0; 
		float diffMouseY = 0; 


		// if a mouse button is pressed allow drag rotation
		if (inputMouseButtons)
		{
			// if the mouse was not pressed last frame. but is pressed this frame, 
			// then we save the old mouse position, button state, and camera rotation
			// to achieve drag rotation
			if (!mouseButtonDownPrevFrame)
			{
				oldMouseX = newMouseX;
				oldMouseY = newMouseY;
				oldAzimuthalAngle = azimuthalAngle;
				oldPolarAngle = polarAngle;
				mouseButtonDownPrevFrame = true;
			}
			else
			{
				// mouse sensitivity
				float mouseXSensitivity = 0.3f;
				float mouseYSensitivity = 0.3f;
				// find the difference between where the mouse is at currently and where the user started
				// dragging it
				diffMouseX = mouseXSensitivity*static_cast<float>(newMouseX - oldMouseX);
				diffMouseY = mouseYSensitivity*static_cast<float>(newMouseY - oldMouseY);
				
				// update rotation based on the mouse difference
				azimuthalAngle = oldAzimuthalAngle + 2.0f*pi*diffMouseX;
				polarAngle = oldPolarAngle + pi * diffMouseY;
			}
		}
		else 
			mouseButtonDownPrevFrame = false;

		// pick forward direction based on the mouse difference from the previous frame
		camera.e2 = vec3(sinf(azimuthalAngle)*cosf(polarAngle), sinf(polarAngle), cosf(azimuthalAngle)*cosf(polarAngle));
		// build an orthonormal system from this forward direction and an up vector
		camera.e0 = normalize(cross(vec3(0, 1, 0), camera.e2));
		camera.e1 = cross(camera.e2, camera.e0);

		// increment frames counter in order to measure fps
		++frames;
	}
	tigrFree(screen);

	//image.savePPM("out.ppm");

	return 0;
}