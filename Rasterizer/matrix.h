#pragma once

#include <iostream>
#include <vector>
#include "vec4.h"
#include "opti.h"

// Matrix class for 4x4 transformation matrices
class matrix {
    union {
        float m[4][4]; // 2D array representation of the matrix
        float a[16];   // 1D array representation of the matrix for linear access
    };

public:
    // Default constructor initializes the matrix as an identity matrix
    matrix() {
        identity();
    }

    // Access matrix elements by row and column
    float& operator()(unsigned int row, unsigned int col) { return m[row][col]; }

    // Display the matrix elements in a readable format
    void display() {
        for (unsigned int i = 0; i < 4; i++) {
            for (unsigned int j = 0; j < 4; j++)
                std::cout << m[i][j] << '\t';
            std::cout << std::endl;
        }
    }

    // Multiply the matrix by a 4D vector
    // Input Variables:
    // - v: vec4 object to multiply with the matrix
    // Returns the resulting transformed vec4
    vec4 operator * (const vec4& v) const {
        vec4 result;
        result[0] = a[0] * v[0] + a[1] * v[1] + a[2] * v[2] + a[3] * v[3];
        result[1] = a[4] * v[0] + a[5] * v[1] + a[6] * v[2] + a[7] * v[3];
        result[2] = a[8] * v[0] + a[9] * v[1] + a[10] * v[2] + a[11] * v[3];
        result[3] = a[12] * v[0] + a[13] * v[1] + a[14] * v[2] + a[15] * v[3];
        return result;
    }

    // Multiply the matrix by another matrix
    // Input Variables:
    // - mx: Another matrix to multiply with
    // Returns the resulting matrix
    matrix operator * (const matrix& mx) const {
        matrix ret;
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                ret.a[row * 4 + col] =
                    a[row * 4 + 0] * mx.a[0 * 4 + col] +
                    a[row * 4 + 1] * mx.a[1 * 4 + col] +
                    a[row * 4 + 2] * mx.a[2 * 4 + col] +
                    a[row * 4 + 3] * mx.a[3 * 4 + col];
            }
        }
        return ret;
    }

    // Create a perspective projection matrix
    // Input Variables:
    // - fov: Field of view in radians
    // - aspect: Aspect ratio of the viewport
    // - n: Near clipping plane
    // - f: Far clipping plane
    // Returns the perspective matrix
    static matrix makePerspective(float fov, float aspect, float n, float f) {
        if (matrixOpti) { // Do some divisions first to make the function as fast as possible.
            matrix m;
            m.zero();
            float tanHalfFov = std::tan(fov / 2.0f);
            float aT = 1.f / (aspect * tanHalfFov);
            float zNorm = 1.f / (f - n);

            m.a[0] = aT;
            m.a[5] = aT * aspect;
            m.a[10] = -f * zNorm;
            m.a[11] = -(f * n) * zNorm;
            m.a[14] = -1.0f;
            return m;
        }
        else {
            matrix m;
            m.zero();
            float tanHalfFov = std::tan(fov / 2.0f);

            m.a[0] = 1.0f / (aspect * tanHalfFov);
            m.a[5] = 1.0f / tanHalfFov;
            m.a[10] = -f / (f - n);
            m.a[11] = -(f * n) / (f - n);
            m.a[14] = -1.0f;
            return m;
        }
    }

    // Create a translation matrix
    // Input Variables:
    // - tx, ty, tz: Translation amounts along the X, Y, and Z axes
    // Returns the translation matrix
    static matrix makeTranslation(float tx, float ty, float tz) {
        matrix m;
        if (!matrixOpti) m.identity(); // This isn't necessary as our default constructor is already identity.
        m.a[3] = tx;
        m.a[7] = ty;
        m.a[11] = tz;
        return m;
    }

    // Create a rotation matrix around the Z-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateZ(float aRad) {
        matrix m;
        if (!matrixOpti) {
            m.identity(); // Unnecessary as default constructor is identity.
            m.a[0] = std::cos(aRad);
            m.a[1] = -std::sin(aRad);
            m.a[4] = std::sin(aRad);
            m.a[5] = std::cos(aRad);
        }
        else { // Calculate cos(aRad) and sin(aRad) first then just reuse it.
            float cosA = std::cos(aRad);
            float sinA = std::sin(aRad);
            m.a[0] = cosA;
            m.a[1] = -sinA;
            m.a[4] = sinA;
            m.a[5] = cosA;
        }
        return m;
    }

    // Create a rotation matrix around the X-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateX(float aRad) {
        matrix m;
        if (!matrixOpti) {
            m.identity(); // Unnecessary as default constructor is identity.
            m.a[5] = std::cos(aRad);
            m.a[6] = -std::sin(aRad);
            m.a[9] = std::sin(aRad);
            m.a[10] = std::cos(aRad);
        }
        else { // Calculate cos(aRad) and sin(aRad) first then just reuse it.
            float cosA = std::cos(aRad);
            float sinA = std::sin(aRad);
            m.a[5] = cosA;
            m.a[6] = -sinA;
            m.a[9] = sinA;
            m.a[10] = cosA;
        }
        return m;
    }

    // Create a rotation matrix around the Y-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateY(float aRad) {
        matrix m;
        if (!matrixOpti) {
            m.identity(); // Unnecessary as default constructor is identity.
            m.a[0] = std::cos(aRad);
            m.a[2] = std::sin(aRad);
            m.a[8] = -std::sin(aRad);
            m.a[10] = std::cos(aRad);
        }
        else { // Calculate cos(aRad) and sin(aRad) first then just reuse it.
            float cosA = std::cos(aRad);
            float sinA = std::sin(aRad);
            m.a[0] = cosA;
            m.a[2] = sinA;
            m.a[8] = -sinA;
            m.a[10] = cosA;
        }
        return m;
    }

    // Create a composite rotation matrix from X, Y, and Z rotations
    // Input Variables:
    // - x, y, z: Rotation angles in radians around each axis
    // Returns the composite rotation matrix
    static matrix makeRotateXYZ(float x, float y, float z) {
        if (matrixOpti) {
            float cx = std::cos(x), sx = std::sin(x);
            float cy = std::cos(y), sy = std::sin(y);
            float cz = std::cos(z), sz = std::sin(z);

            matrix m;

            m.a[0] = cy * cz;
            m.a[1] = -cy * sz;
            m.a[2] = sy;
            m.a[3] = 0.0f;

            m.a[4] = sx * sy * cz + cx * sz;
            m.a[5] = -sx * sy * sz + cx * cz;
            m.a[6] = -sx * cy;
            m.a[7] = 0.0f;

            m.a[8] = -cx * sy * cz + sx * sz;
            m.a[9] = cx * sy * sz + sx * cz;
            m.a[10] = cx * cy;
            m.a[11] = 0.0f;

            m.a[12] = 0.0f;
            m.a[13] = 0.0f;
            m.a[14] = 0.0f;
            m.a[15] = 1.0f;

            return m;
        }

        return matrix::makeRotateX(x) * matrix::makeRotateY(y) * matrix::makeRotateZ(z);
    }

    // Create a scaling matrix
    // Input Variables:
    // - s: Scaling factor
    // Returns the scaling matrix
    static matrix makeScale(float s) {
        matrix m;
        s = std::max(s, 0.01f); // Ensure scaling factor is not too small
        if (!matrixOpti) m.identity(); // Unnecessary as default constructor is identity.
        m.a[0] = s;
        m.a[5] = s;
        m.a[10] = s;
        return m;
    }

    // Create an identity matrix
    // Returns an identity matrix
    static matrix makeIdentity() {
        matrix m;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m.m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
        return m;
    }

private:
    // Set all elements of the matrix to 0
    void zero() {
        for (unsigned int i = 0; i < 16; i++)
            a[i] = 0.f;
    }

    // Set the matrix as an identity matrix
    void identity() {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }
};


