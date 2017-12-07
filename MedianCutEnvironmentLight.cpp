
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



#include "stdafx.h"
#include "lights/MedianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"
//#define MID_DEBUG	
//#define SEG_IMG_DEBUG

//int LightCount=0;
// MedianCutEnvironmentLight Utility Classes
struct MedianCutEnvironmentCube {
    // InfiniteAreaCube Public Methods
    MedianCutEnvironmentCube(const MedianCutEnvironmentLight *l, const Scene *s,
							float t, bool cv, float pe)
							: light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};



// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
	//int width = 0, height =0;
	//LightNum = 32;
	//CutLight = new PointLight[256];

    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
    delete[] texels;
    // Initialize sampling PDFs for infinite area light


	illumination = new float* [height];
	red = new float* [height];
	green = new float* [height];
	blue = new float* [height];

	for(int i=0;i<height;i++){
		illumination[i] = new float [width];
		red[i] = new float [width];
		green[i] = new float [width];
		blue[i] = new float [width];
	}

    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
			//Accumulate the illumination
			if(v==0 && u==0){
				illumination[0][0] = img[0];
				red[0][0] = radianceMap->Lookup(up, vp, filter).CallColor(0) * sinTheta;
				green[0][0] = radianceMap->Lookup(up, vp, filter).CallColor(1) * sinTheta;
				blue[0][0] = radianceMap->Lookup(up, vp, filter).CallColor(2) * sinTheta;
			}
			else if(v==0 && u>0){
				illumination[v][u] = illumination[v][u-1] + img[u];
				red[v][u] = illumination[v][u-1] + radianceMap->Lookup(up, vp, filter).CallColor(0) * sinTheta;
				green[v][u] = illumination[v][u-1] + radianceMap->Lookup(up, vp, filter).CallColor(1) * sinTheta;
				blue[v][u] = illumination[v][u-1] + radianceMap->Lookup(up, vp, filter).CallColor(2) * sinTheta;
			}
			else if(u==0 && v>0){
				illumination[v][u] = illumination[v-1][u] + img[v*width];
				red[v][u] = illumination[v-1][u] + radianceMap->Lookup(up, vp, filter).CallColor(0) * sinTheta;
				green[v][u] = illumination[v-1][u] + radianceMap->Lookup(up, vp, filter).CallColor(1) * sinTheta;
				blue[v][u] = illumination[v-1][u] + radianceMap->Lookup(up, vp, filter).CallColor(2) * sinTheta;
			}
			else{
				illumination[v][u] = illumination[v-1][u] + illumination[v][u-1] + img[u+v*width] - illumination[v-1][u-1];
				red[v][u] = illumination[v-1][u] + illumination[v][u-1] + radianceMap->Lookup(up, vp, filter).CallColor(0) * sinTheta - illumination[v-1][u-1];
				green[v][u] = illumination[v-1][u] + illumination[v][u-1] + radianceMap->Lookup(up, vp, filter).CallColor(1) * sinTheta - illumination[v-1][u-1];
				blue[v][u] = illumination[v-1][u] + illumination[v][u-1] + radianceMap->Lookup(up, vp, filter).CallColor(2) * sinTheta - illumination[v-1][u-1];
			}
			//if(v==height-1) printf("%10.2lf",illumination[v][u]);		
        }
    }
	//system("pause");
	//Seperate the longest edge into two section by illumination
	//determine the longest edge
	Segment temp = { illumination[height-1][width-1],{0,0,0}, width, height, 0, 0, 0, 0, false };
	ImageSegment.push_back(temp);


	for( int i=0; i<nSamples; i++ ) {
	if( ImageSegment[ImageSegment.size()-1].width >= ImageSegment[ImageSegment.size()-1].height )
			Dir = Horizon;
		else
			Dir = Vertical;


		//Find seperating position-----------
		//printf("\nHalfIllumination Start\n");
		float HalfIllumination = SubImgIllumination( ImageSegment[ImageSegment.size()-1].left, ImageSegment[ImageSegment.size()-1].left + ImageSegment[ImageSegment.size()-1].width-1,
										ImageSegment[ImageSegment.size()-1].top, ImageSegment[ImageSegment.size()-1].top + ImageSegment[ImageSegment.size()-1].height-1,
										width, height, illumination );
		//printf("HalfIllumination Finish\n");
		//printf("Origin:\nleft:%d width:%d top:%d height:%d Illumination:%10.2lf ImgW:%d ImgH:%d\n",ImageSegment[ImageSegment.size()-1].left,ImageSegment[ImageSegment.size()-1].width,ImageSegment[ImageSegment.size()-1].top,ImageSegment[ImageSegment.size()-1].height,HalfIllumination,width, height);
		HalfIllumination /= 2;
		int mid = FindIlluminationPosition( ImageSegment[ImageSegment.size()-1].left, ImageSegment[ImageSegment.size()-1].left + ImageSegment[ImageSegment.size()-1].width-1,
											ImageSegment[ImageSegment.size()-1].top, ImageSegment[ImageSegment.size()-1].top + ImageSegment[ImageSegment.size()-1].height-1,
											width, height, Dir, illumination, HalfIllumination );

		if(  (Dir==Horizon  && (mid==ImageSegment[ImageSegment.size()-1].left || mid==(ImageSegment[ImageSegment.size()-1].left + ImageSegment[ImageSegment.size()-1].width-1) )) ||
			 (Dir==Vertical && (mid==ImageSegment[ImageSegment.size()-1].top  || mid==(ImageSegment[ImageSegment.size()-1].top + ImageSegment[ImageSegment.size()-1].height-1) ))
		  ){
			CannotSegment.push_back(ImageSegment[ImageSegment.size()-1]);
			ImageSegment.pop_back();
			i--;
		}
		else{
		//Insert new segments and delete the used one
			InsertSegment( &ImageSegment, mid, width, height, Dir, illumination, red, green, blue );	
		}
	}
	for(int i=0;i<CannotSegment.size();i++)
		ImageSegment.push_back(CannotSegment[i]);

	// Compute sampling distributions for rows and columns of image
	//float *SegImg;
	//for(int i=0;i<nSamples;i++){
	//	SegImg = new float [ ImageSegment[i].height*ImageSegment[i].width ];
	//	for(int row=0;row<ImageSegment[i].height;row++){
	//		for(int col=0;col<ImageSegment[i].width;col++)
	//			SegImg[col+row*width] = img[ (col+ImageSegment[i].left)+(row+ImageSegment[i].top)*width];
	//	}
	//	distribution.push_back( new Distribution2D(SegImg, ImageSegment[i].width, ImageSegment[i].height) );
	//	delete[] SegImg;
	//}

	distribution = new Distribution2D(img, width, height);
	// Create definite Point Light 
	GiveAllSegmentCenter( &ImageSegment, width, height, illumination );

	CannotSegment.clear();
    delete[] img;
	for(int i=0;i<height;i++){
		delete [] illumination[i];
		delete [] red[i];
		delete [] green[i];
		delete [] blue[i];
	}
	delete [] illumination;
	delete [] red;
	delete [] green;
	delete [] blue;
	//ImageSegment.clear();
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _MedianCutEnvironmentLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _MedianCutEnvironmentLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _MedianCutEnvironmentLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _MedianCutEnvironmentLight_ to SH from cube map sampling
        SHProjectCube(MedianCutEnvironmentCube(this, scene, time, computeLightVis,
												pEpsilon),  
					  p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_MEDIAN_CUT_ENV_LIGHT_STARTED_SAMPLE();
    float uv[2];
	int randomLightIndex = (int)ls.uPos[0]*(nSamples-1);
	uv[0] = ImageSegment[randomLightIndex].c_x/width;
	uv[1] = ImageSegment[randomLightIndex].c_y/height;

    // Convert infinite light sample point to direction
    float theta =  uv[1]* M_PI, phi =  uv[0]* 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));
    // Compute PDF for sampled infinite light direction
    *pdf = 1./(float)nSamples;

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
	Spectrum Ls = Spectrum(ImageSegment[randomLightIndex].color);

	PBRT_MEDIAN_CUT_ENV_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_MEDIAN_CUT_ENV_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) / 
           (2.f * M_PI * M_PI * sintheta);
    PBRT_MEDIAN_CUT_ENV_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_MEDIAN_CUT_ENV_LIGHT_STARTED_SAMPLE();
	//float uv[2],mapPdf;
	//int randomLightIndex = (int)ls.uPos[0]*(nSamples-1);
	//Spectrum Ls;
	////if( LightCount<nSamples ){ 
	//uv[0] = ImageSegment[randomLightIndex].c_x/width;
	//uv[1] = ImageSegment[randomLightIndex].c_y/height;

	//distribution->SampleContinuous(uv[0], uv[1], uv, &mapPdf);

	//float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
 //   float costheta = cosf(theta), sintheta = sinf(theta);
 //   float sinphi = sinf(phi), cosphi = cosf(phi);
 //   Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
 //                                   costheta));
 //   *Ns = (Normal)d;

 //   // Compute origin for infinite light sample ray
 //   Point worldCenter;
 //   float worldRadius;
 //   scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
 //   Vector v1, v2;
 //   CoordinateSystem(-d, &v1, &v2);
 //   float d1, d2;
 //   ConcentricSampleDisk(u1, u2, &d1, &d2);
 //   Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
 //   *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);
 //
 //   // Compute _MedianCutEnvironmentLight_ ray PDF
 //   float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
 //   float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
 //   *pdf = directionPdf * areaPdf;
 //   if (sintheta == 0.f) *pdf = 0.f;
 //   Ls = (radianceMap->Lookup(ImageSegment[randomLightIndex].c_x, ImageSegment[randomLightIndex].c_y), SPECTRUM_ILLUMINANT);
 //   
	//printf("scene randomLightIndex:%d\n",randomLightIndex);
	////LightCount += 1;
	////}
	////else
	//	//return Spectrum(0);
   // Compute direction for infinite light sample ray

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	// Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _InfiniteAreaLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
	PBRT_MEDIAN_CUT_ENV_LIGHT_FINISHED_SAMPLE();
    return Ls;
}




double MedianCutEnvironmentLight::SubImgIllumination( int left, int right, int top, int bottom, int &ImageWidth, int &ImageHeight, float** illumination){
	//try{
		if( left<= right && top<=bottom && bottom<ImageHeight && right<ImageWidth && left>=0 && top>=0 ){
			if( left==0 && top!=0 )
				return ( illumination[bottom][right] - illumination[top-1][right]);
			else if( top==0 && left!=0 )
				return  ( illumination[bottom][right]- illumination[bottom][left-1]);
			else if( top==0 && left==0 )
				return  illumination[bottom][right];
			else
				return  ( illumination[bottom][right] - illumination[top-1][right] - illumination[bottom][left-1] + illumination[top-1][left-1] );
		}
		else{
			//throw NewException("There\'s wrong input \'SubImgIllumination\' ");
		printf("There\'s wrong input \'SubImgIllumination\' \n");
		printf("left:%d right:%d top:%d bottom:%d ImgW:%d, ImgH:%d\n",left,right,top,bottom,ImageWidth,ImageHeight);
		system("pause");
		exit(1);
		}
	//}
	//catch(NewException e){
		//std::cerr<< e.message()<< std::endl;
		//printf("There\'s wrong input \'SubImgIllumination\' \n");
		//printf("left:%d right:%d top:%d bottom:%d ImgW:%d, ImgH:%d\n",left,right,top,bottom,ImageWidth,ImageHeight);
		//system("pause");
		//exit(1);
	//}
}

//Find the seperating position 
int MedianCutEnvironmentLight::FindIlluminationPosition( int &left, int right, int &top, int bottom, int &ImageWidth, int &ImageHeight,IlluminationDirection dir, float** illumination,float TargetIllumination){
	double MidIllumination;
	int low, high, mid;
	//binary search
	switch(dir){
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	case Horizon:
		low = left;  high = right;  mid = 0;
		while (low <= high)
		{
			mid = (low + high) / 2;
			if ( (top==0 && left==0 && illumination[bottom][mid] == TargetIllumination) ||
				 (top==0 && left!=0 && (illumination[bottom][mid]-illumination[bottom][left-1]) == TargetIllumination) ||
				 (top!=0 && left==0 && (illumination[bottom][mid]-illumination[top-1][mid]) == TargetIllumination) ||
				 (top!=0 && left!=0 && (illumination[bottom][mid]-illumination[top-1][mid]-illumination[bottom][left-1] + illumination[top-1][left-1]) == TargetIllumination) 
				)
			{
#ifdef MID_DEBUG
			//double MidIllumination = (illumination[ImageSegment[ImageSegment.size()-1].top+ImageSegment[ImageSegment.size()-1].height-1][mid]-illumination[ImageSegment[ImageSegment.size()-1].top-1][mid]-illumination[ImageSegment[ImageSegment.size()-1].top+ImageSegment[ImageSegment.size()-1].height-1][ImageSegment[ImageSegment.size()-1].left-1] + illumination[ImageSegment[ImageSegment.size()-1].top-1][ImageSegment[ImageSegment.size()-1].left-1]);
			printf("Insert Hor:\nmid:%d Illumination:%10.2lf\n",mid, TargetIllumination );		
			system("pause");
#endif
				return mid;
			}
			else if( (top==0 && left==0 && illumination[bottom][mid] > TargetIllumination) ||
					 (top==0 && left!=0 && (illumination[bottom][mid]-illumination[bottom][left-1]) > TargetIllumination) ||
					 (top!=0 && left==0 && (illumination[bottom][mid]-illumination[top-1][mid]) > TargetIllumination) ||
					 (top!=0 && left!=0 && (illumination[bottom][mid]-illumination[top-1][mid]-illumination[bottom][left-1] + illumination[top-1][left-1]) > TargetIllumination) 
					)
			{
				high = mid - 1;
			}
			else if( (top==0 && left==0 && illumination[bottom][mid] < TargetIllumination) ||
					 (top==0 && left!=0 && (illumination[bottom][mid]-illumination[bottom][left-1]) < TargetIllumination) ||
					 (top!=0 && left==0 && (illumination[bottom][mid]-illumination[top-1][mid]) < TargetIllumination) ||
					 (top!=0 && left!=0 && (illumination[bottom][mid]-illumination[top-1][mid]-illumination[bottom][left-1] + illumination[top-1][left-1]) < TargetIllumination) 
					)
			{
				low = mid + 1;
			}
			else;
		}
#ifdef MID_DEBUG
			if(top==0 && left==0) MidIllumination = illumination[bottom][mid];
			else if(top==0 && left!=0 ) MidIllumination = (illumination[bottom][mid]-illumination[bottom][left-1]);
			else if(top!=0 && left==0 ) MidIllumination = (illumination[bottom][mid]-illumination[top-1][mid]);
			else if(top!=0 && left!=0 ) MidIllumination = (illumination[bottom][mid]-illumination[top-1][mid]-illumination[bottom][left-1] + illumination[top-1][left-1]);
			printf("Insert Hor:\nmid:%d Illumination:%10.2lf\n",mid, MidIllumination );		
			system("pause");
#endif
		return mid;
		break;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	case Vertical:
		low = top;  high = bottom;  mid = 0;
		while (low <= high)
		{
			mid = (low + high) / 2;
			if ( (left==0 && top==0 && illumination[mid][right] == TargetIllumination) ||
				 (left==0 && top!=0 && (illumination[mid][right]-illumination[top-1][right]) == TargetIllumination) ||
				 (left!=0 && top==0 && (illumination[mid][right]-illumination[mid][left-1]) == TargetIllumination) ||
				 (left!=0 && top!=0 && (illumination[mid][right]-illumination[mid][left-1]-illumination[top-1][right]+illumination[top-1][left-1]) == TargetIllumination) 
				)
			{
#ifdef MID_DEBUG
			//double MidIllumination = (illumination[ImageSegment[ImageSegment.size()-1].top+ImageSegment[ImageSegment.size()-1].height-1][mid]-illumination[ImageSegment[ImageSegment.size()-1].top-1][mid]-illumination[ImageSegment[ImageSegment.size()-1].top+ImageSegment[ImageSegment.size()-1].height-1][ImageSegment[ImageSegment.size()-1].left-1] + illumination[ImageSegment[ImageSegment.size()-1].top-1][ImageSegment[ImageSegment.size()-1].left-1]);
			printf("Insert Ver:\nmid:%d Illumination:%10.2lf\n",mid, TargetIllumination );		
			system("pause");
#endif
				return mid;
			}
			else if ( (left==0 && top==0 && illumination[mid][right] > TargetIllumination) ||
					  (left==0 && top!=0 && (illumination[mid][right]-illumination[top-1][right]) > TargetIllumination) ||
					  (left!=0 && top==0 && (illumination[mid][right]-illumination[mid][left-1]) > TargetIllumination) ||
					  (left!=0 && top!=0 && (illumination[mid][right]-illumination[mid][left-1]-illumination[top-1][right]+illumination[top-1][left-1]) > TargetIllumination) 
					 )
			{
				high = mid - 1;
			}
			else if ( (left==0 && top==0 && illumination[mid][right] < TargetIllumination) ||
					  (left==0 && top!=0 && (illumination[mid][right]-illumination[top-1][right]) < TargetIllumination) ||
					  (left!=0 && top==0 && (illumination[mid][right]-illumination[mid][left-1]) < TargetIllumination) ||
					  (left!=0 && top!=0 && (illumination[mid][right]-illumination[mid][left-1]-illumination[top-1][right]+illumination[top-1][left-1]) < TargetIllumination) 
					)
			{
				low = mid + 1;
			}
			else;
		}
#ifdef MID_DEBUG
			if(left==0 && top==0 ) MidIllumination = illumination[mid][right];
			else if(left==0 && top!=0 ) MidIllumination = (illumination[mid][right]-illumination[top-1][right]);
			else if(left!=0 && top==0 ) MidIllumination = (illumination[mid][right]-illumination[mid][left-1]);
			else if(left!=0 && top!=0 ) MidIllumination = (illumination[mid][right]-illumination[mid][left-1]-illumination[top-1][right]+illumination[top-1][left-1]);
			printf("Insert Ver:\nmid:%d Illumination:%10.2lf\n",mid, MidIllumination );		
			system("pause");
#endif
		return mid;
		break;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	default:
		break;
	}
}

void MedianCutEnvironmentLight::InsertSegment( vector<Segment> *Img, int &mid, int &ImageWidth, int &ImageHeight, IlluminationDirection dir, float **illumination, float **red, float **green, float **blue ){
	//Set two SubSegments--------------------
#ifdef SEG_IMG_DEBUG
	printf("origin: left:%d width:%d top:%d height:%d illu:%10.2lf\n",(*Img)[(*Img).size()-1].left,(*Img)[(*Img).size()-1].width,(*Img)[(*Img).size()-1].top,(*Img)[(*Img).size()-1].height,(*Img)[(*Img).size()-1].illumination);
	(dir==0)? ( printf("Horizon\n") ):( printf("Vertical\n") );
	printf("\n");
	if( (*Img).size()==38 ){
		printf("mid:%d\n",mid);
	}
#endif	
	Segment temp[2];
	switch(dir){
	case Horizon:
		//left
		temp[0].left = (*Img)[(*Img).size()-1].left;						temp[0].top = (*Img)[(*Img).size()-1].top;
		temp[0].width = mid+1-temp[0].left;									temp[0].height = (*Img)[(*Img).size()-1].height ;
		temp[0].c_x = 0;													temp[0].c_y = 0;
		temp[0].IsCenterRight = false;
		//printf("Insert Hor1 Start\n");
		temp[0].illumination = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, illumination );
		temp[0].color[0] = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, red );
		temp[0].color[1] = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, green );
		temp[0].color[2] = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, blue );
		//printf("Insert Hor1 Finish\n");
		//right
		temp[1].left = mid+1;												temp[1].top = (*Img)[(*Img).size()-1].top;
		temp[1].width = (*Img)[(*Img).size()-1].width - temp[0].width;		temp[1].height = (*Img)[(*Img).size()-1].height ;
		temp[1].c_x = 0;													temp[1].c_y = 0;
		temp[1].IsCenterRight = false;
		//printf("Insert Hor2 Start\n");
		temp[1].illumination = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, illumination );
		temp[1].color[0] = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, red );
		temp[1].color[1] = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, green );
		temp[1].color[2] = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, blue );
		//printf("Insert Hor2 Finish\n");
		break;
	case Vertical:
		//upward
		temp[0].left = (*Img)[(*Img).size()-1].left;						temp[0].top = (*Img)[(*Img).size()-1].top;
		temp[0].width =(*Img)[(*Img).size()-1].width;						temp[0].height = mid+1-temp[0].top;
		temp[0].c_x = 0;													temp[0].c_y = 0;
		temp[0].IsCenterRight = false;
		//printf("Insert Ver1 Start\n");
		temp[0].illumination = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, illumination );
		temp[0].color[0] = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, red );
		temp[0].color[1] = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, green );
		temp[0].color[2] = SubImgIllumination( temp[0].left, temp[0].left + temp[0].width-1, temp[0].top, temp[0].top + temp[0].height-1, ImageWidth, ImageHeight, blue );
		
		//printf("Insert Ver1 Finish\n");
		//down
		temp[1].left = (*Img)[(*Img).size()-1].left;						temp[1].top = temp[0].top + temp[0].height;
		temp[1].width = (*Img)[(*Img).size()-1].width;						temp[1].height = (*Img)[(*Img).size()-1].height - temp[0].height;
		temp[1].c_x = 0;													temp[1].c_y = 0;
		temp[1].IsCenterRight = false;
		//printf("Insert Ver2 Start\n");
		temp[1].illumination = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, illumination );
		temp[1].color[0] = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, red );
		temp[1].color[1] = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, green );
		temp[1].color[2] = SubImgIllumination( temp[1].left, temp[1].left + temp[1].width-1, temp[1].top, temp[1].top + temp[1].height-1, ImageWidth, ImageHeight, blue );
		//printf("Insert Ver2 Finish\n");
		break;
	default:
		break;
	}

	//Insert---------------------------------
	(*Img).pop_back(); //Clear the first element
	if((*Img).empty()){
		if( temp[0].illumination <= temp[1].illumination ){
			(*Img).insert( (*Img).begin(),temp[1] );
			(*Img).insert( (*Img).begin(),temp[0] );}
		else{
			(*Img).insert( (*Img).begin(),temp[0] );
			(*Img).insert( (*Img).begin(),temp[1] );}
	}
	else{
			int index;
			bool IsInsert;
			for(int i=0; i<2; i++){
				index = 0;
				IsInsert = false;
				while( index<(*Img).size() ){
					if( temp[i].illumination > (*Img)[index].illumination ){
						index++;
						//printf("i:%d Index:%d\n",i,index);
					}
					else{
						(*Img).insert( (*Img).begin()+index , temp[i] );
						IsInsert = true;
						break;
					}
				}
				if(!IsInsert) (*Img).insert( (*Img).begin()+index , temp[i] );
			}
			//int index = InsertIndexBinarySearch( 0, (*Img).size()-1, Img, temp[0].illumination );		
			//if( (*Img)[index].illumination>=temp[0].illumination )			//If the value of the index is 'bigger' than INSERTED value 
			//	(*Img).insert( (*Img).begin()+(index),temp[0] );			//Insert at index																	
			//else															//If the value of the index is 'smaller' than INSERTED value 
			//	(*Img).insert( (*Img).begin()+(index+1),temp[0] );			//Insert at index+1
			//index = InsertIndexBinarySearch( 0, (*Img).size()-1, Img, temp[1].illumination );
			//if( (*Img)[index].illumination>=temp[1].illumination )			//If the value of the index is 'bigger' than INSERTED value 
			//	(*Img).insert( (*Img).begin()+(index),temp[1] );			//Insert at index
			//else															//If the value of the index is 'smaller' than INSERTED value 
			//	(*Img).insert( (*Img).begin()+(index+1),temp[1] );			//Insert at index+1
			
	}
#ifdef SEG_IMG_DEBUG
	for(int i=0;i<(*Img).size();i++)
		printf("%d: left:%d width:%d top:%d height:%d illu:%10.2lf\n",i,(*Img)[i].left,(*Img)[i].width,(*Img)[i].top,(*Img)[i].height,(*Img)[i].illumination);
	printf("\n");
	system("pause");
#endif

}


bool MedianCutEnvironmentLight::GiveAllSegmentCenter(vector<Segment> *Img, int &ImageWidth, int &ImageHeight, float** illumination){
	double InvertIllumination, PixelIllumination;
	double temp_cx, temp_cy;
	bool AllCenterIsRight = true;

	for(int i=0;i<(*Img).size();i++){
		InvertIllumination = 1/(*Img)[i].illumination;
		temp_cx = temp_cy = 0;
		for(int row=0; row<(*Img)[i].height; row++ ){
			for(int col=0; col<(*Img)[i].width; col++ ){
				PixelIllumination = SubImgIllumination( (*Img)[i].left + col, (*Img)[i].left + col, (*Img)[i].top + row, (*Img)[i].top + row, ImageWidth, ImageHeight, illumination );
				temp_cx += ( PixelIllumination*( (*Img)[i].left + col )*InvertIllumination );
				temp_cy += ( PixelIllumination*( (*Img)[i].top + row  )*InvertIllumination );
			}
		}
		(*Img)[i].c_x = temp_cx / ImageWidth;
		(*Img)[i].c_y = temp_cy / ImageHeight;
		(*Img)[i].IsCenterRight = true;
		AllCenterIsRight &= (*Img)[i].IsCenterRight;
	}

	return AllCenterIsRight;
}
