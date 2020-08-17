#include "pch.h"
#include <iostream>
#include <algorithm>
#include <windows.h>
#include <omp.h>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "our_gl.h"

Model *model = NULL;
const int width = 800;
const int height = 800;
const int depth = 255;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);

Vec3f light_dir(1, 1, 1);
Vec3f       eye(1, 1, 3);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);

struct GouraudShader : public IShader {
	Vec3f varying_intensity;
	virtual Vec4f vertex(int iface, int nthvert) {
		varying_intensity[nthvert] = (std::max)(0.0f, model->normal(iface, nthvert)*light_dir);
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));
		return Viewport * Projection * ModelView * gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color)
	{
		float intensity = varying_intensity * bar;
		color = TGAColor(255, 255, 255)*intensity;
		return false;
	}
};
//void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(p0.x - p1.x) < std::abs(p0.y- p1.y))
//	{
//		std::swap(p0.x, p0.y);
//		std::swap(p1.x, p1.y);
//		steep = true;
//	}
//	if (p0.x > p1.x)
//	{
//		std::swap(p0.x, p1.x);
//		std::swap(p0.y, p1.y);
//	}
//	int dx = p1.x - p0.x;
//	int dy = p1.y - p0.y;
//	float derror2 = std::abs(dy)*2;
//	float error2 = 0;
//	int y = p0.y;
//	for (int x = p0.x; x <= p1.x; x++)
//	{
//		if (steep)
//		{
//			image.set(y, x, color);
//		}
//		else
//		{
//			image.set(x, y, color);
//		}
//		error2 += derror2;
//		if (error2 > dx)
//		{
//			y += (p1.y > p0.y ? 1 : -1);
//			error2 -= dx*2;
//		}
//	}
//}
//扫描线填充
//void triangle(Vec2i t0,Vec2i t1,Vec2i t2,TGAImage &image,TGAColor color)
//{
//	if (t0.y == t1.y && t0.y == t2.y) return;
//	if (t0.y > t1.y) std::swap(t0, t1);
//	if (t0.y > t2.y) std::swap(t0, t2);
//	if (t1.y > t2.y) std::swap(t1, t2);
//	int total_height = t2.y - t0.y;
//	for (int i = 0; i <total_height; i++)
//	{
//		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
//		int segment_height = second_half ? t2.y-t1.y:t1.y-t0.y;
//		float alpha = (float) i / total_height;
//		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here 
//		Vec2i A = t0 + (t2 - t0)*alpha;
//		Vec2i B = second_half ? t1 + (t2 - t1)*beta : t0 + (t1 - t0)*beta;
//		if (A.x > B.x) std::swap(A, B);
//		for (int j = A.x; j <= B.x; j++) {
//			image.set(j, t0.y+i, color); // attention, due to int casts t0.y+i != A.y 
//		}
//	}
//}


//edge equation
//bool edge_equation(Vec2i P0, Vec2i P1, Vec2i P2, Vec2i P)
//{
//	int eq_cur = (P2.y - P1.y)*P.x + (P1.x - P2.x)*P.y + (P2.x*P1.y - P1.x * P2.y);
//	int eq_p0 = (P2.y - P1.y)*P0.x + (P1.x - P2.x)*P0.y + (P2.x*P1.y - P1.x * P2.y);
//	if (eq_cur >= 0 && eq_p0 >= 0)return true;
//	if (eq_cur <= 0 && eq_p0 <= 0) return true;
//	return false;
//}

//void rasterize(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color, int ybuffer[])
//{
//	if (p0.x > p1.x)
//	{
//		std::swap(p0, p1);
//	}
//	for (int x = p0.x; x < p1.x; x++)
//	{
//		float t = (x - p0.x) / (float)(p1.x - p0.x);
//		int y = p0.y*(1.0 - t) + p1.y*t;
//		if (ybuffer[x] < y)
//		{
//			ybuffer[x] = y;
//			for(int i=0;i<image.get_height();i++)
//				image.set(x, i, color);
//		}
//	}
//}

int main()
{
	DWORD dwTimeStarted;
	// model import
	model = new Model("obj/african_head.obj");	

	lookat(eye, center, up);
	viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
	projection(-1.f / (eye - center).norm());
	light_dir.normalize();
	
	TGAImage image(width, height, TGAImage::RGB);
	TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

	dwTimeStarted = ::GetTickCount();
    #pragma omp parallel num_threads(6)
	GouraudShader shader;
	for (int i = 0; i < model->nfaces(); i++)
	{
		Vec4f screen_coords[3];
		for (int j = 0; j < 3; j++)
		{
			screen_coords[j] = shader.vertex(i, j);
		}
		triangle(screen_coords, shader, image, zbuffer);
	}
	image.flip_vertically(); // to place the origin in the bottom left corner of the image
	zbuffer.flip_vertically();
	image.write_tga_file("output-shader.tga");
	zbuffer.write_tga_file("zbuffer-shader.tga");

	delete model;
	return 0;
}


