#pragma once

#include <cmath>
#include <cstdint>
#include "opti.h"

// The `colour` class represents an RGB colour with floating-point precision.
// It provides various utilities for manipulating and converting colours.
class colour { // r,g,b,a?
    union {
        struct {
            float r, g, b; // Red, Green, and Blue components of the colour
        };
        float rgb[3];     // Array representation of the RGB components
    };

public:
    // Enum for indexing the RGB components
    enum Colour { RED = 0, GREEN = 1, BLUE = 2 };

    // Constructor to initialize the colour with specified RGB values.
    // Default values are 0 (black).
    // Input Variables:
    // - _r: Red component (default 0.0f)
    // - _g: Green component (default 0.0f)
    // - _b: Blue component (default 0.0f)
    colour(float _r = 0, float _g = 0, float _b = 0) : r(_r), g(_g), b(_b) {}

    // Sets the RGB components of the colour.
    // Input Variables:
    // - _r: Red component
    // - _g: Green component
    // - _b: Blue component
    void set(float _r, float _g, float _b) { r = _r, g = _g, b = _b; }

    // Accesses the specified component of the colour by index.
    // Input Variables:
    // - c: Index of the component (RED, GREEN, or BLUE)
    // Returns a reference to the specified component.
    float& operator[] (Colour c) { return rgb[c]; }

    // Assigns the values of another colour to this one.
    // Input Variables:
    // - c: The source color
    void operator = (colour c) { // CONST COLOUR C?
        r = c.r;
        g = c.g;
        b = c.b;
    }

    // Clamps the RGB components of the colour to the range [0, 1].
    void clampColour() {
        if (colourOpti) { // Same thing as bottom but without using std::min
            if (r > 1.0f) r = 1.0f; else if (r < 0.0f) r = 0.0f;
            if (g > 1.0f) g = 1.0f; else if (g < 0.0f) g = 0.0f;
            if (b > 1.0f) b = 1.0f; else if (b < 0.0f) b = 0.0f;
        }
        else {
            r = std::min(r, 1.0f);
            g = std::min(g, 1.0f);
            b = std::min(b, 1.0f);
        }
    }

    // Converts the floating-point RGB values to integer values (0-255).
    // Output Variables:
    // - cr: Red component as an unsigned char
    // - cg: Green component as an unsigned char
    // - cb: Blue component as an unsigned char
    void toRGB(unsigned char& cr, unsigned char& cg, unsigned char& cb) {
        if (colourOpti) { // Cast to int instead of floor. There's no difference for positive numbers. If we clamp before calling toRGB its the same functionality.
            cr = static_cast<unsigned char>((int)(r * 255));
            cg = static_cast<unsigned char>((int)(g * 255));
            cb = static_cast<unsigned char>((int)(b * 255));
        } else {
            cr = static_cast<unsigned char>(std::floor(r * 255));
            cg = static_cast<unsigned char>(std::floor(g * 255));
            cb = static_cast<unsigned char>(std::floor(b * 255));
        }
    }

    // Scales the RGB components of the colour by a scalar value.
    // Input Variables:
    // - scalar: The scaling factor
    // Returns a new `colour` object with scaled components.
    colour operator * (const float scalar) {
        if (colourOpti) return colour(r * scalar, g * scalar, b * scalar); // No need to create a new object in memory.
        colour c;
        c.r = r * scalar;
        c.g = g * scalar;
        c.b = b * scalar;
        return c;
    }

    // Multiplies the RGB components of this colour with another colour.
    // Input Variables:
    // - col: The other color to multiply with
    // Returns a new `colour` object with multiplied components.
    colour operator * (const colour& col) {
        if (colourOpti) return colour(r * col.r, g * col.g, b * col.b); // No need to create a new object in memory.
        colour c;
        c.r = r * col.r;
        c.g = g * col.g;
        c.b = b * col.b;
        return c;
    }

    // Adds the RGB components of another colour to this one.
    // Input Variables:
    // - _c: The other colour to add
    // Returns a new `colour` object with added components.
    colour operator + (const colour& _c) {
        if (colourOpti) return colour(r + _c.r, g + _c.g, b + _c.b); // No need to create a new object in memory.
        colour c;
        c.r = r + _c.r;
        c.g = g + _c.g;
        c.b = b + _c.b;
        return c;
    }
};
