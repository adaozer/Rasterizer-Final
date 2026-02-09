#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include "opti.h"

// Simple support class for a 2D vector
class vec2D {
public:
    float x, y;

    // Default constructor initializes both components to 0
    vec2D() { x = y = 0.f; };

    // Constructor initializes components with given values
    vec2D(float _x, float _y) : x(_x), y(_y) {}

    // Constructor initializes components from a vec4
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }

    // Display the vector components
    void display() { std::cout << x << '\t' << y << std::endl; }

    // Overloaded subtraction operator for vector subtraction
    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// Class representing a triangle for rendering purposes
class triangle {
    Vertex v[3];       // Vertices of the triangle
    float area;        // Area of the triangle
    colour col[3];     // Colors for each vertex of the triangle

public:
    // Constructor initializes the triangle with three vertices
    // Input Variables:
    // - v1, v2, v3: Vertices defining the triangle
    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        // Calculate the 2D area of the triangle
        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = std::fabs(e1.x * e2.y - e1.y * e2.x);
    }

    // Helper function to compute the cross product for barycentric coordinates
    // Input Variables:
    // - v1, v2: Edges defining the vector
    // - p: Point for which coordinates are being calculated
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    // Compute barycentric coordinates for a given point
    // Input Variables:
    // - p: Point to check within the triangle
    // Output Variables:
    // - alpha, beta, gamma: Barycentric coordinates of the point
    // Returns true if the point is inside the triangle, false otherwise
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma, float a = 1) {
        if (triangleOpti) {
            alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) * a; // Multiply by 1/area instead of dividing every time
            if (alpha < 0.f) return false; // Check after we get each component so we can stop early if necessary.
            beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) * a;
            if (beta < 0.f) return false;
            gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) * a;
            if (gamma < 0.f) return false;
        }
        else {
            alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
            beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
            gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

            if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        }

        return true;
    }

    // Template function to interpolate values using barycentric coordinates
    // Input Variables:
    // - alpha, beta, gamma: Barycentric coordinates
    // - a1, a2, a3: Values to interpolate
    // Returns the interpolated value
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    // Draw the triangle on the canvas
    // Input Variables:
    // - renderer: Renderer object for drawing
    // - L: Light object for shading calculations
    // - ka, kd: Ambient and diffuse lighting coefficients
    void draw(Renderer& renderer, Light& L, float ka, float kd, int clipY_start = 0, int clipY_end = INT32_MAX) {
        if (triangleOpti) {
            float signedArea = (v[1].p[0] - v[0].p[0]) * (v[2].p[1] - v[0].p[1]) - (v[1].p[1] - v[0].p[1]) * (v[2].p[0] - v[0].p[0]);
            if (signedArea <= 0.0f) return;

            int canvasW = renderer.canvas.getWidth();
            int canvasH = renderer.canvas.getHeight();

            float minX = std::min({ v[0].p[0], v[1].p[0], v[2].p[0] });
            float maxX = std::max({ v[0].p[0], v[1].p[0], v[2].p[0] });
            float minY = std::min({ v[0].p[1], v[1].p[1], v[2].p[1] });
            float maxY = std::max({ v[0].p[1], v[1].p[1], v[2].p[1] });

            int x0 = std::max(0, std::min(canvasW - 1, (int)std::floor(minX)));
            int x1 = std::max(0, std::min(canvasW - 1, (int)std::ceil(maxX)));
            int y0 = std::max(0, std::max(clipY_start, std::min(canvasH - 1, (int)std::floor(minY))));
            int y1 = std::max(0, std::min(clipY_end, std::min(canvasH - 1, (int)std::ceil(maxY))));

            if (x0 >= x1 || y0 >= y1) return;

            vec4 omega = L.omega_i;
            omega.normalise();
            colour ambientTerm = L.ambient * ka;
            colour lightIntensity = L.L;

            vec2D A(v[0].p), B(v[1].p), C(v[2].p);

            vec2D pRow((float)x0 + 0.5f, (float)y0 + 0.5f);

            float E0_row = getC(A, B, pRow);
            float E1_row = getC(B, C, pRow);
            float E2_row = getC(C, A, pRow);

            vec2D e0 = B - A; vec2D e1 = C - B; vec2D e2 = A - C;
            float E0_dx = -e0.y, E0_dy = e0.x;
            float E1_dx = -e1.y, E1_dy = e1.x;
            float E2_dx = -e2.y, E2_dy = e2.x;
            float invArea = 1.f / area;
            for (int y = y0; y < y1; ++y) {
                float E0 = E0_row, E1 = E1_row, E2 = E2_row;
                for (int x = x0; x < x1; ++x) {

                    if (E0 >= 0.f && E1 >= 0.f && E2 >= 0.f) {
                        float alpha = E0 * invArea;
                        float beta = E1 * invArea;
                        float gamma = E2 * invArea;

                        float depth = (v[0].p[2] * beta) + (v[1].p[2] * gamma) + (v[2].p[2] * alpha);
                        float& zBuff = renderer.zbuffer(x, y);

                        if (depth > 0.001f && zBuff > depth) {
                            zBuff = depth;

                            vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
                            normal.normalise();

                            float dot = std::max(vec4::dot(omega, normal), 0.0f);
                            colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);

                            colour out = (c * kd) * (lightIntensity * dot) + ambientTerm;

                            unsigned char r, g, b;
                            out.clampColour();
                            out.toRGB(r, g, b);
                            renderer.canvas.draw(x, y, r, g, b);
                        }
                    }
                    E0 += E0_dx; E1 += E1_dx; E2 += E2_dx;
                }
                E0_row += E0_dy; E1_row += E1_dy; E2_row += E2_dy;
            }
        }
        else {
            vec2D minV, maxV;

            // Get the screen-space bounds of the triangle
            getBoundsWindow(renderer.canvas, minV, maxV);

            // Skip very small triangles
            if (area < 1.f) return;
            float invArea = 1.f;

            // Iterate over the bounding box and check each pixel
            for (int y = (int)(minV.y); y < (int)ceil(maxV.y); y++) {
                for (int x = (int)(minV.x); x < (int)ceil(maxV.x); x++) {
                    float alpha, beta, gamma;

                    // Check if the pixel lies inside the triangle
                    if (getCoordinates(vec2D((float)x, (float)y), alpha, beta, gamma)) {
                        // Interpolate color, depth, and normals
                        colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
                        c.clampColour();
                        float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
                        vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
                        normal.normalise();

                        // Perform Z-buffer test and apply shading
                        if (renderer.zbuffer(x, y) > depth && depth > 0.001f) {
                            // typical shader begin
                            L.omega_i.normalise();
                            float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                            colour a = (c * kd) * (L.L * dot) + (L.ambient * ka); // using kd instead of ka for ambient
                            // typical shader end
                            unsigned char r, g, b;
                            a.toRGB(r, g, b);
                            renderer.canvas.draw(x, y, r, g, b);
                            renderer.zbuffer(x, y) = depth;
                        }
                    }
                }
            }
        }
    }
    // Compute the 2D bounds of the triangle
    // Output Variables:
    // - minV, maxV: Minimum and maximum bounds in 2D space
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = std::min(minV.x, v[i].p[0]);
            minV.y = std::min(minV.y, v[i].p[1]);
            maxV.x = std::max(maxV.x, v[i].p[0]);
            maxV.y = std::max(maxV.y, v[i].p[1]);
        }
    }

    // Compute the 2D bounds of the triangle, clipped to the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    // Output Variables:
    // - minV, maxV: Clipped minimum and maximum bounds
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = std::max(minV.x, static_cast<float>(0));
        minV.y = std::max(minV.y, static_cast<float>(0));
        maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
        maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
    }

    // Debugging utility to display the triangle bounds on the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};
