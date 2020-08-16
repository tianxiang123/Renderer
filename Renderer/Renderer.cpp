#include "pch.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <windows.h>
#include <omp.h>
#include "tgaimage.h"
#include "model.h"

const int width = 800;
const int height = 800;
const int depth = 255;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);

Vec3f light_dir(0, 0, -1);
Vec3f camera(0, 0, 3);

Vec3f m2v(Matrix m) {
	return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
}

Matrix v2m(Vec3f v) {
	Matrix m(4, 1);
	m[0][0] = v.x;
	m[1][0] = v.y;
	m[2][0] = v.z;
	m[3][0] = 1.f;
	return m;
}

Matrix viewport(int x, int y, int w, int h) {
	Matrix m = Matrix::identity(4);
	m[0][3] = x + w / 2.f;
	m[1][3] = y + h / 2.f;
	m[2][3] = depth / 2.f;

	m[0][0] = w / 2.f;
	m[1][1] = h / 2.f;
	m[2][2] = depth / 2.f;
	return m;
}

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
	bool steep = false;
	if (std::abs(p0.x - p1.x) < std::abs(p0.y- p1.y))
	{
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x > p1.x)
	{
		std::swap(p0.x, p1.x);
		std::swap(p0.y, p1.y);
	}
	int dx = p1.x - p0.x;
	int dy = p1.y - p0.y;
	float derror2 = std::abs(dy)*2;
	float error2 = 0;
	int y = p0.y;
	for (int x = p0.x; x <= p1.x; x++)
	{
		if (steep)
		{
			image.set(y, x, color);
		}
		else
		{
			image.set(x, y, color);
		}
		error2 += derror2;
		if (error2 > dx)
		{
			y += (p1.y > p0.y ? 1 : -1);
			error2 -= dx*2;
		}
	}
}
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

//重心坐标计算
Vec3f barycentric(Vec3f *pts, Vec3f P)
{
	Vec3f u = Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]) ^ Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]);
	if (std::abs(u[2]) > 1e-2 ) return Vec3f(1.0f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vec3f(-1, 1, 1);
}

//edge equation
//bool edge_equation(Vec2i P0, Vec2i P1, Vec2i P2, Vec2i P)
//{
//	int eq_cur = (P2.y - P1.y)*P.x + (P1.x - P2.x)*P.y + (P2.x*P1.y - P1.x * P2.y);
//	int eq_p0 = (P2.y - P1.y)*P0.x + (P1.x - P2.x)*P0.y + (P2.x*P1.y - P1.x * P2.y);
//	if (eq_cur >= 0 && eq_p0 >= 0)return true;
//	if (eq_cur <= 0 && eq_p0 <= 0) return true;
//	return false;
//}
void triangle(Vec3f *pts, Vec2f* texts,float *zbuffer, TGAImage &image, Model* model) {
	Vec2f bboxmin((std::numeric_limits<float>::max)(), (std::numeric_limits<float>::max)());
	Vec2f bboxmax(-(std::numeric_limits<float>::max)(), -(std::numeric_limits<float>::max)());
	Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			bboxmin[j] = (std::max)(0.0f, (std::min)(bboxmin[j], pts[i][j]));
			bboxmax[j] = (std::min)(clamp[j], (std::max)(bboxmax[j], pts[i][j]));
		}
	}
	Vec3f P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
	{
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
		{
			//重心坐标
			Vec3f bc_screen = barycentric(pts, P); 
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			P.z = 0;
			Vec2f Ptexcoord(0, 0);
			for (int i = 0; i < 3; i++)
			{
				P.z += pts[i][2] * bc_screen[i];
				Ptexcoord[0] += texts[i][0] * bc_screen[i];
				Ptexcoord[1] += texts[i][1] * bc_screen[i];
			}
				
			if (zbuffer[int(P.x + P.y*width)] < P.z)
			{
				zbuffer[int(P.x + P.y*width)] = P.z;
				TGAColor color = model->diffuse(Ptexcoord);
				image.set(P.x, P.y, color);
			}
			//edge equotion
			//if (edge_equation(pts[0], pts[1], pts[2], P) && edge_equation(pts[1], pts[2], pts[0], P) &&
			//	edge_equation(pts[2], pts[0], pts[1], P))
			//{
			//	image.set(P.x, P.y, color);
			//}
		}
	}
}

void rasterize(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color, int ybuffer[])
{
	if (p0.x > p1.x)
	{
		std::swap(p0, p1);
	}
	for (int x = p0.x; x < p1.x; x++)
	{
		float t = (x - p0.x) / (float)(p1.x - p0.x);
		int y = p0.y*(1.0 - t) + p1.y*t;
		if (ybuffer[x] < y)
		{
			ybuffer[x] = y;
			for(int i=0;i<image.get_height();i++)
				image.set(x, i, color);
		}
	}
}

Vec3f world2screen(Vec3f v) {
	return Vec3f(int((v.x + 1.)*width / 2. + .5), int((v.y + 1.)*height / 2. + .5), v.z);
}

int main()
{
	DWORD dwTimeStarted;
	TGAImage diffuseimg(width, height, TGAImage::RGB);
	// model import
	Model* model = new Model("obj/african_head.obj");	

	float *zbuffer = new float[width*height];
	for (int i = width * height; i--; zbuffer[i] = -(std::numeric_limits<float>::max)());

	dwTimeStarted = ::GetTickCount();
#pragma omp parallel num_threads(12)
	{ // draw the model
		Matrix Projection = Matrix::identity(4);
		Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
		Projection[3][2] = -1.f / camera.z;

		TGAImage image(width, height, TGAImage::RGB);
		for (int i = 0; i < model->nfaces(); i++) {
			std::vector<int> face = model->face(i);
			Vec3f screen_coords[3];
			Vec3f world_coords[3];
			for (int j = 0; j < 3; j++) {
				Vec3f v = model->vert(face[j]);
				screen_coords[j] = m2v(ViewPort*Projection*v2m(v));
				world_coords[j] = v;
			}
			Vec3f n = (world_coords[2] - world_coords[0])^(world_coords[1] - world_coords[0]);
			n.normalize();
			float intensity = n * light_dir;
			if (intensity > 0) {
				Vec2f uv[3];
				for (int k = 0; k < 3; k++) {
					uv[k] = model->uv(i, k);
				}
				triangle(screen_coords, uv, zbuffer,image,model);
			}
		}

		image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
	}

	{ // dump z-buffer (debugging purposes only)
		TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				zbimage.set(i, j, TGAColor(zbuffer[i + j * width]*255, zbuffer[i + j * width] * 255, zbuffer[i + j * width] * 255,255));
			}
		}
		zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		zbimage.write_tga_file("zbuffer.tga");
	}
	//for (int i = 0; i < model->nfaces(); i++) {
	//	std::vector<int> face = model->face(i);
	//	Vec2i screen_coords[3];
	//	Vec3f world_coords[3];
	//	for (int j = 0; j < 3; j++) {
	//		Vec3f v = model->vert(face[j]);
	//		screen_coords[j] = Vec2i((v.x+1.0)*width/2.0,(v.y+1.0)*height/2.0);
	//		world_coords[j] = v;
	//	}
	//	Vec3f n =cross(world_coords[2] - world_coords[0],world_coords[1] - world_coords[0]);
	//	n.normalize();
	//	float intensity = n * lightDir;
	//	if (intensity > 0)
	//	{
	//		triangle(screen_coords, frame, TGAColor(intensity * 255, intensity * 255, intensity * 255));
	//	}	
	//}
	//for (int i = 0; i < model->nfaces(); i++)
	//{
	//	std::vector<int> face = model->face(i);
	//	Vec3f pts[3];
	//	Vec3f world_coords[3];
	//	Vec2f texts[3];
	//	for (int j = 0; j < 3; j++) {
	//		world_coords[j] = model->vert(face[j]);
	//		pts[j] = world2screen(model->vert(face[j]));
	//		texts[j] = model->uv(i, j);
	//	}
	//	Vec3f n =cross(world_coords[2] - world_coords[0],world_coords[1] - world_coords[0]);
	//	n.normalize();
	//	float intensity = n * lightDir;
	//	triangle(pts, texts,zbuffer, diffuseimg, model);
	//}
	
	


	diffuseimg.flip_vertically();
	diffuseimg.write_tga_file("zBuffer.tga");
	delete zbuffer;
	delete model;
	std::cout << ::GetTickCount() - dwTimeStarted;
	return 0;
}


