/*#include<iostream>
#include<pcl/io/pcd_io.h>
#include<pcl/point_types.h>
#include<pcl/point_cloud.h>
#include<pcl/kdtree/kdtree_flann.h>
#include<pcl/kdtree/io.h>
#include<pcl/search/kdtree.h>
#include <pcl/search/search.h>
#include <pcl/features/normal_3d.h>
#include <pcl/common/common.h>
#include <pcl/common/pca.h>
#include<vector>
#include<algorithm>
#include<fstream>
#include<limits.h>
#include "fftw3.h"
#include <numbers>
#include<cstdlib>
#include<cmath>
#include <iomanip>
#include <math.h>
#pragma comment(lib, "libfftw3-3.lib")
using namespace std;
#define REAL 0
#define IMAG 1
std::vector<float> g_kix(8);
std::vector<float> g_kiy(8);

float k_i(vector<float> SampV, int N)
{
    //cout << "input your sample value......" << endl;
    //vector<double> SampValue; //採樣值

    vector<float> a; //FFT後實部運算
    vector<float> b; //FFT後虛部運算
    fftw_complex* x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); //fftw專用輸出資料結構(輸入)
    fftw_complex* y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); //(輸出)
    for (int i = 0; i < N; i++) {
        x[i][REAL] = SampV[i]; //採樣值皆令其在實部
        x[i][IMAG] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, x, y, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    for (int i = 0; i < N; i++) {
        if (i == 0) {
            a.push_back(y[i][REAL] / N);
            b.push_back(0);
        }
        else {
            a.push_back((2 * y[i][REAL]) / N);
            b.push_back((2 * y[i][IMAG]) / N);
        }
    }

    float var_1 = 0;
    float var_2 = 0;
    for (int n = 0; n < N; n++) {
        var_1 = var_1 + pow(n + 1, 2) * a[n];
        var_2 = var_2 + (n + 1) * b[n];
    }

    float result = -(pow(M_PI, 2) / pow(2 * M_PI, 2) * var_1 / pow((1 + pow(var_2 / 2, 2)), 1.5));

    fftw_destroy_plan(plan);

    return result;
}

void kix(int index, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr newcloud, pcl::KdTreeFLANN<pcl::PointXYZ> newkdtree, Eigen::Matrix3d vector, int* ptr) {


    for (int u = 1; u < 9; u++) {
        int q1index = -1;
        int q2index = -1;
        int q3index = -1;
        int q4index = -1;

        pcl::PointXYZ newsearchPoint;

        newsearchPoint.x = (*cloud)[index].x + u * (vector(ptr[0])) / 100;
        newsearchPoint.y = (*cloud)[index].y + u * (vector(ptr[0] + 1)) / 100;
        newsearchPoint.z = (*cloud)[index].z + u * (vector(ptr[0] + 2)) / 100;

        int newK = 30;

        std::vector<int> newpointIdxNKNSearch(newK);
        std::vector<float> newpointNKNSquaredDistance(newK);

        float XPZa = vector(ptr[1]) * newsearchPoint.x;
        float XPZb = vector(ptr[1] + 1) * newsearchPoint.y;
        float XPZc = vector(ptr[1] + 2) * newsearchPoint.z;

        float XPYa = vector(ptr[2]) * newsearchPoint.x;
        float XPYb = vector(ptr[2] + 1) * newsearchPoint.y;
        float XPYc = vector(ptr[2] + 2) * newsearchPoint.z;

        if (newkdtree.nearestKSearch(newsearchPoint, newK, newpointIdxNKNSearch, newpointNKNSquaredDistance) > 0)
        {


            for (std::size_t i = 0; i < newpointIdxNKNSearch.size(); ++i) {

                float XOZa = vector(ptr[1]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float XOZb = vector(ptr[1] + 1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float XOZc = vector(ptr[1] + 2) * (*newcloud)[newpointIdxNKNSearch[i]].z;

                float XOYa = vector(ptr[2]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float XOYb = vector(ptr[2] + 1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float XOYc = vector(ptr[2] + 2) * (*newcloud)[newpointIdxNKNSearch[i]].z;

                if (XOZa + XOZb + XOZc < XPZa + XPZb + XPZc && XOYa + XOYb + XOYc < XPYa + XPYb + XPYc && q1index < 0) {
                    q1index = newpointIdxNKNSearch[i];

                }
                else if (XOZa + XOZb + XOZc < XPZa + XPZb + XPZc && XOYa + XOYb + XOYc > XPYa + XPYb + XPYc && q2index < 0) {
                    q2index = newpointIdxNKNSearch[i];
                }
                else if (XOZa + XOZb + XOZc > XPZa + XPZb + XPZc && XOYa + XOYb + XOYc < XPYa + XPYb + XPYc && q3index < 0) {
                    q3index = newpointIdxNKNSearch[i];
                }
                else if (XOZa + XOZb + XOZc > XPZa + XPZb + XPZc && XOYa + XOYb + XOYc > XPYa + XPYb + XPYc && q4index < 0) {
                    q4index = newpointIdxNKNSearch[i];
                }
            }

        }

        std::vector<int> qindex;
        qindex.push_back(q1index);
        qindex.push_back(q2index);
        qindex.push_back(q3index);
        qindex.push_back(q4index);

        float u1down = 0;
        float g_ui = 0;

        for (int j = 0; j < 4; j++)
        {
            if (qindex.at(j) >= 0) {
                u1down += 1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(j)].x, 2) +
                    pow(newsearchPoint.y - newcloud->points[qindex.at(j)].y, 2) +
                    pow(newsearchPoint.z - newcloud->points[qindex.at(j)].z, 2)));
            }
        }
        if (ptr[2] == 0) {
            for (int i = 0; i < 4; i++)
            {
                if (qindex.at(i) >= 0) {
                    g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                        pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                        pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].x) / u1down;
                }
            }
        }
        else if (ptr[2] == 3)
        {
            for (int i = 0; i < 4; i++)
            {
                if (qindex.at(i) >= 0) {
                    g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                        pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                        pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].y) / u1down;
                }
            }
        }
        else if (ptr[2] == 6)
        {
            for (int i = 0; i < 4; i++)
            {
                if (qindex.at(i) >= 0) {
                    g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                        pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                        pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].z) / u1down;
                }
            }
        }

        g_kix[u - 1] = g_ui;

    }
}

void kiy(int index, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr newcloud, pcl::KdTreeFLANN<pcl::PointXYZ> newkdtree, Eigen::Matrix3d vector, int* ptr) {

    for (int u = 1; u < 9; u++) {
        int q1index = -1;
        int q2index = -1;
        int q3index = -1;
        int q4index = -1;


        pcl::PointXYZ newsearchPoint;

        newsearchPoint.x = (*cloud)[index].x + u * (vector(ptr[1])) / 100;
        newsearchPoint.y = (*cloud)[index].y + u * (vector(ptr[1] + 1)) / 100;
        newsearchPoint.z = (*cloud)[index].z + u * (vector(ptr[1] + 2)) / 100;

        int newK = 30;

        std::vector<int> newpointIdxNKNSearch(newK);
        std::vector<float> newpointNKNSquaredDistance(newK);

        float YPZa = vector(ptr[0]) * newsearchPoint.x;
        float YPZb = vector(ptr[0] + 1) * newsearchPoint.y;
        float YPZc = vector(ptr[0] + 2) * newsearchPoint.z;

        float XPYa = vector(ptr[2]) * newsearchPoint.x;
        float XPYb = vector(ptr[2] + 1) * newsearchPoint.y;
        float XPYc = vector(ptr[2] + 2) * newsearchPoint.z;


        if (newkdtree.nearestKSearch(newsearchPoint, newK, newpointIdxNKNSearch, newpointNKNSquaredDistance) > 0)
        {

            for (std::size_t i = 0; i < newpointIdxNKNSearch.size(); ++i) {

                float YOZa = vector(ptr[0]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float YOZb = vector(ptr[0] + 1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float YOZc = vector(ptr[0] + 2) * (*newcloud)[newpointIdxNKNSearch[i]].z;

                float XOYa = vector(ptr[2]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float XOYb = vector(ptr[2] + 1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float XOYc = vector(ptr[2] + 2) * (*newcloud)[newpointIdxNKNSearch[i]].z;




                if (YOZa + YOZb + YOZc < YPZa + YPZb + YPZc && XOYa + XOYb + XOYc < XPYa + XPYb + XPYc && q1index < 0) {
                    q1index = newpointIdxNKNSearch[i];

                }
                else if (YOZa + YOZb + YOZc < YPZa + YPZb + YPZc && XOYa + XOYb + XOYc > XPYa + XPYb + XPYc && q2index < 0) {
                    q2index = newpointIdxNKNSearch[i];
                }
                else if (YOZa + YOZb + YOZc > YPZa + YPZb + YPZc && XOYa + XOYb + XOYc < XPYa + XPYb + XPYc && q3index < 0) {
                    q3index = newpointIdxNKNSearch[i];
                }
                else if (YOZa + YOZb + YOZc > YPZa + YPZb + YPZc && XOYa + XOYb + XOYc > XPYa + XPYb + XPYc && q4index < 0) {
                    q4index = newpointIdxNKNSearch[i];
                }
            }

        }




        std::vector<int> qindex;
        qindex.push_back(q1index);
        qindex.push_back(q2index);
        qindex.push_back(q3index);
        qindex.push_back(q4index);

        float u1down = 0;
        float g_ui = 0;

        for (int j = 0; j < 4; j++)
        {
            if (qindex.at(j) >= 0) {
                u1down += 1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(j)].x, 2) +
                    pow(newsearchPoint.y - newcloud->points[qindex.at(j)].y, 2) +
                    pow(newsearchPoint.z - newcloud->points[qindex.at(j)].z, 2)));
            }
        }
        if (ptr[2] == 0) {
            for (int i = 0; i < 4; i++)
            {
                if (qindex.at(i) >= 0) {
                    g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                        pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                        pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].x) / u1down;
                }
            }
        }
        else if (ptr[2] == 3)
        {
            for (int i = 0; i < 4; i++)
            {
                if (qindex.at(i) >= 0) {
                    g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                        pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                        pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].y) / u1down;
                }
            }
        }
        else if (ptr[2] == 6)
        {
            for (int i = 0; i < 4; i++)
            {
                if (qindex.at(i) >= 0) {
                    g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                        pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                        pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].z) / u1down;
                }
            }
        }

        g_kiy[u - 1] = g_ui;

    }
}*/