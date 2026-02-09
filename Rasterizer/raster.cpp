#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"
#include "opti.h"
#include <mutex>
#include <thread>
#include <atomic>
#include <condition_variable>

// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.

struct BakedTri { // Includes all the data required for rasterisation. No need to deal with meshes/matrices etc.
    Vertex v[3];
    float ka, kd; 
};

class RenderSystem {
public:
    static RenderSystem& getInstance() {
        static RenderSystem instance; // Singleton
        return instance;
    }

    // This is the manager so threads aren't created again each frame.
    Renderer* currentRenderer = nullptr;
    Light* currentLight = nullptr;
    std::vector<BakedTri> drawList; // List of triangles in the scene

    // For threading
    std::vector<std::thread> workers;
    std::mutex mtx;
    std::condition_variable cv_worker;
    std::condition_variable cv_main;

    bool shutdown = false;
    int currentFrame = 0;
    std::vector<int> threadFrame;
    std::atomic<int> threadsFinished{ 0 };
    int numThreads;


    RenderSystem() {
        drawList.reserve(20000);
        numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;

        threadFrame.resize(numThreads, 0);

        for (int i = 0; i < numThreads; i++) {
            workers.emplace_back(&RenderSystem::workerLoop, this, i);
        }
    }

    ~RenderSystem() {
        {
            std::unique_lock<std::mutex> lock(mtx);
            shutdown = true;
        }

        cv_worker.notify_all();
        for (auto& t : workers) t.join();
    }

    void workerLoop(int id) {
        while (true) {
            std::unique_lock<std::mutex> lock(mtx);

            cv_worker.wait(lock, [this, id] {
                return currentFrame > threadFrame[id] || shutdown;
                });

            if (shutdown) return;

            threadFrame[id] = currentFrame; // if new frame,
            lock.unlock();

            int height = currentRenderer->canvas.getHeight(); // Find your slice of the screen
            int sliceHeight = height / numThreads;
            int startY = id * sliceHeight;
            int endY = (id == numThreads - 1) ? height : (id + 1) * sliceHeight;

            for (const auto& triData : drawList) {  // Draw triangles within their slice of the screen

                float minY = std::min({ triData.v[0].p[1], triData.v[1].p[1], triData.v[2].p[1] });
                float maxY = std::max({ triData.v[0].p[1], triData.v[1].p[1], triData.v[2].p[1] });

                if (maxY < startY || minY > endY) continue;

                triangle t(triData.v[0], triData.v[1], triData.v[2]);
                t.draw(*currentRenderer, *currentLight, triData.ka, triData.kd, startY, endY);
            }

            int finished = ++threadsFinished;
            if (finished == numThreads) {
                cv_main.notify_one();
            }
        }
    }

    void processFrame(Renderer& renderer, Light& L) {
        currentRenderer = &renderer; // move to next frame if the current frame is done
        currentLight = &L;
        threadsFinished = 0;

        {
            std::unique_lock<std::mutex> lock(mtx);
            currentFrame++;
        }
        cv_worker.notify_all();

        {
            std::unique_lock<std::mutex> lock(mtx);
            cv_main.wait(lock, [this] { return threadsFinished == numThreads; });
        }
    }
};

void renderScene(Renderer& renderer, std::vector<Mesh*>& scene, matrix& camera, Light& L) {
    // preprocessing
    // this prepares the scene for the thread workers
    RenderSystem& sys = RenderSystem::getInstance();
    sys.drawList.clear();

    float halfW = renderer.canvas.getWidth() * 0.5f;
    float halfH = renderer.canvas.getHeight() * 0.5f;

    static std::vector<Vertex> tempVerts;

    for (Mesh* mesh : scene) {
        matrix mvp = renderer.perspective * camera * mesh->world;
        matrix world = mesh->world;

        if (tempVerts.size() < mesh->vertices.size()) {
            tempVerts.resize(mesh->vertices.size());
        }

        size_t vCount = mesh->vertices.size();
        for (size_t i = 0; i < vCount; ++i) {
            Vertex& outV = tempVerts[i];
            const Vertex& inV = mesh->vertices[i];
            outV.p = mvp * inV.p;
            outV.normal = world * inV.normal;
            outV.normal.normalise();
            outV.rgb = inV.rgb;
        }

        for (const triIndices& ind : mesh->triangles) {
            if (tempVerts[ind.v[0]].p[3] < 0.1f ||
                tempVerts[ind.v[1]].p[3] < 0.1f ||
                tempVerts[ind.v[2]].p[3] < 0.1f) continue;

            BakedTri bTri;
            bTri.ka = mesh->ka;
            bTri.kd = mesh->kd;

            for (int k = 0; k < 3; k++) {
                bTri.v[k] = tempVerts[ind.v[k]];

                bTri.v[k].p.divideW();

                bTri.v[k].p[0] = (bTri.v[k].p[0] * halfW) + halfW;
                bTri.v[k].p[1] = halfH - (bTri.v[k].p[1] * halfH);
            }

            float area = (bTri.v[1].p[0] - bTri.v[0].p[0]) * (bTri.v[2].p[1] - bTri.v[0].p[1]) -
                (bTri.v[1].p[1] - bTri.v[0].p[1]) * (bTri.v[2].p[0] - bTri.v[0].p[0]);

            if (area <= 0) continue;

            sys.drawList.push_back(bTri);
        }
    }

    sys.processFrame(renderer, L);
}

void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;

    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal; 
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
            t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        // Create a triangle object and render it
        triangle tri(t[0], t[1], t[2]);
        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix::makeIdentity();
    }
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() // custom scene, basically copy/mesh of scene 1 and scene 2.
{
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    bool running = true;
    std::vector<Mesh*> scene;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    struct rRot { float x; float y; float z; };
    std::vector<rRot> rotations;

    const unsigned int rings = 18;
    const float ringStepZ = 3.5f;
    const float ringRadius = 3.0f;
    const unsigned int cubesPerRing = 8;

    for (unsigned int i = 0; i < rings; i++) {
        float z = -(ringStepZ * (float)i);

        for (unsigned int k = 0; k < cubesPerRing; k++) {
            float a = (2.0f * (float)M_PI) * ((float)k / (float)cubesPerRing);
            float x = std::cos(a) * ringRadius;
            float y = std::sin(a) * ringRadius;

            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            m->world = matrix::makeTranslation(x, y, z) * makeRandomRotation();
            scene.push_back(m);

            rRot r{ rng.getRandomFloat(-.08f, .08f),
                    rng.getRandomFloat(-.08f, .08f),
                    rng.getRandomFloat(-.08f, .08f) };
            rotations.push_back(r);
        }
    }

    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);

    float sphereOffset = -4.0f;
    float sphereStep = 0.08f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -8.f);

    float zoffset = 6.0f;
    float step = -0.08f;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        camera = matrix::makeTranslation(0, 0, -zoffset);
        zoffset += step;
        if (zoffset < -55.f || zoffset > 6.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        for (unsigned int i = 0; i < rotations.size(); i++) {
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
        }

        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -8.f);
        if (sphereOffset > 4.0f || sphereOffset < -4.0f) {
            sphereStep *= -1.f;
        }

        renderScene(renderer, scene, camera, L);
        renderer.present();
    }

    for (auto& m : scene) delete m;
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    bool running = true;

    std::vector<Mesh*> scene;

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (multiThread) {
            renderScene(renderer, scene, camera, L);
        }
        else {
            for (auto& m : scene)
                render(renderer, m, camera, L);
        }
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    std::vector<Mesh*> scene;

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.push_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        if (multiThread) {
            renderScene(renderer, scene, camera, L); // For multithreading, we use a different render function that takes the scene as param rather than individual meshes.
        }
        else {
            for (auto& m : scene)
                render(renderer, m, camera, L);
        }
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Entry point of the application
// No input variables
int main() {
    // Uncomment the desired scene function to run
    //scene1();
    //scene2();
    sceneTest(); 
    

    return 0;
}