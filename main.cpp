#include <iostream>
#include "vec.h"
#include "sphere.h"
#include "camera.h"
#include "scene.h"
#include "integrator.h"
#include "light.h"

using namespace HydrogenCG;
int main()
{
	Image image;
	image.init(1024, 768);
	Camera camera;

	// Please use modern C++, std::unique_ptr should be used here.
	Scene scene;
	scene.surfaces.push_back(new Sphere({ 0, 0, 3 }, 1, { 1, 0.3f, 0 }));
	scene.surfaces.push_back(new Sphere({ -1.5f, 0, 5 }, 1, { 0.2f, 0.7f, 0.1f }));
	scene.surfaces.push_back(new Sphere({ 0, -1001, 3 }, 1000, vec3(0.4f)));
	scene.lights.push_back(LightPoint({ 0, 3, 0 }, vec3(30.0f)));

	// I think integral math is over the head of most beginners.
	// Imo we should use a term here that's still fitting the concept,
	// but beginners will actually understand.
	// Also, please don't leak memory. Beginners will likely use our tutorials
	// as the foundations of their projects, and we should not teach them bad habits.
	// I propose the use of "auto" here as follows: auto integrator = IntegratorDiffuseLocal;
	Integrator* integrator = new IntegratorDiffuseLocal;
	integrator->render(image, camera, scene);

	image.savePPM("out.ppm");

	return 0;
}
