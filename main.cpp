#include <iostream>
#include "vec.hpp"
#include "sphere.hpp"
#include "camera.hpp"
#include "scene.hpp"
#include "integrator.hpp"
#include "light.hpp"


using namespace HydrogenCG;
int main()
{
	Image image;
	image.init(640, 480);
	Camera camera;
	Scene scene;
	scene.surfaces.push_back(new Sphere({ 0, 0, 3 }, 1, { 1, 0.3f, 0 }));
	scene.surfaces.push_back(new Sphere({ -1.5f, 0, 5 }, 1, { 0.2f, 0.7f, 0.1f }));
	scene.surfaces.push_back(new Sphere({ 0, -1001, 3 }, 1000, vec3(0.2f, 0.4f, 0.8f)));
	scene.lights.push_back(LightPoint({ 0, 3, 0 }, vec3(10.0f)));
	const Integrator& integrator = IntegratorDiffuse{};
	integrator.render(image, camera, scene);

	image.savePPM("out.ppm");

	return 0;
}