
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_MEDIAN_CUT_ENVIRONMENT_LIGHT_H
#define PBRT_LIGHTS_MEDIAN_CUT_ENVIRONMENT_LIGHT_H

// lights/infinite.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"
//#include "point.h"

enum IlluminationDirection{ Horizon, Vertical };
struct Segment{
	float illumination,color[3];//r g b
	int width, height;
	int left, top;
	float c_x,c_y;
	bool IsCenterRight;
};

class NewException{
public:
	NewException(const char* message) : _message(message){}
	std::string message(){ return _message; }
protected:
	std::string _message;
};

// MedianCutEnvironmentLight Declarations
class MedianCutEnvironmentLight : public Light {
public:
    // MedianCutEnvironmentLight Public Methods
    MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
        const string &texmap);
    ~MedianCutEnvironmentLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;

	double SubImgIllumination( int left, int right, int top, int bottom, int& ImageWidth, int& ImageHeight, float** illumination);
	int FindIlluminationPosition( int& left, int right, int& top, int bottom, int& ImageWidth, int& ImageHeight,IlluminationDirection dir, float** illumination,float TargetIllumination);
	void InsertSegment( vector<Segment> *Img, int& mid, int &ImageWidth, int &ImageHeight, IlluminationDirection dir, float **illumination, float **red, float **green, float **blue );
	//int InsertIndexBinarySearch( int low, int high,vector<Segment> *Img, double& SegmentIllumination );
	bool GiveAllSegmentCenter(vector<Segment> *Img, int &ImageWidth, int &ImageHeight, float** illumination);
private:
    // MedianCutEnvironmentLight Private Data
    MIPMap<RGBSpectrum> *radianceMap;
	Distribution2D *distribution;
    //vector<Distribution2D *> distribution;

	float **illumination;
	float **red, **green, **blue;
	IlluminationDirection Dir;
	vector<Segment> ImageSegment;
	vector<Segment> CannotSegment;
	//int LightNum;
	int width, height;
};


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
														   const ParamSet &paramSet);

#endif // PBRT_LIGHTS_MEDIAN_CUT_ENVIRONMENT_LIGHT_H
