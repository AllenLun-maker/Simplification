#include<iostream>
#include<pcl/io/pcd_io.h>
#include<pcl/point_types.h>
#include<pcl/point_cloud.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/cloud_viewer.h>
#include<pcl/kdtree/kdtree_flann.h>
#include<pcl/kdtree/io.h>
#include<pcl/search/kdtree.h>
#include <pcl/search/search.h>
#include <pcl/features/normal_3d.h>
#include <pcl/common/common.h>
#include <pcl/common/pca.h>
#include <pcl/features/integral_image_normal.h>  //法線估計類標頭檔案
#include <pcl/io/ply_io.h>
#include <pcl/console/parse.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/passthrough.h>
#include <pcl/segmentation/region_growing.h>
#include<vector>
#include<algorithm>
#include<fstream>
#include<limits.h>
#include "fftw3.h"
#include <numbers>
#include<cstdlib>
#include<cmath>
#include <iomanip>
#define _USE_MATH_DEFINES // for C++
#include <math.h>
#pragma comment(lib, "libfftw3-3.lib")

using namespace std;

#define PI 3.14 //定義PI=3.14159
const float Min = INT_MIN; //常數 表無限大

std::vector<float> g_kix(8);
std::vector<float> g_kiy(8);

typedef pcl::PointXYZ PointT;
typedef pcl::PointXYZRGB PointRGB;
typedef pcl::Normal Normal;
#define REAL 0
#define IMAG 1
//實部與虛部

Eigen::Matrix3d eigen_value(string dataname) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::io::loadPCDFile(dataname + ".pcd", *cloud);
    int cld_sz_1 = cloud->size();
    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    normals->resize(cld_sz_1);

    //计算中心点坐标
    double center_x1 = 0, center_y1 = 0, center_z1 = 0;
    for (int i = 0; i < cld_sz_1; i++) {
        center_x1 += cloud->points[i].x;
        center_y1 += cloud->points[i].y;
        center_z1 += cloud->points[i].z;

    }
    center_x1 /= cld_sz_1;
    center_y1 /= cld_sz_1;
    center_z1 /= cld_sz_1;
    //计算协方差矩阵
    double xx1 = 0, xy1 = 0, xz1 = 0, yy1 = 0, yz1 = 0, zz1 = 0;
    for (int i = 0; i < cld_sz_1; i++) {
        xx1 += (cloud->points[i].x - center_x1) * (cloud->points[i].x - center_x1);
        xy1 += (cloud->points[i].x - center_x1) * (cloud->points[i].y - center_y1);
        xz1 += (cloud->points[i].x - center_x1) * (cloud->points[i].z - center_z1);
        yy1 += (cloud->points[i].y - center_y1) * (cloud->points[i].y - center_y1);
        yz1 += (cloud->points[i].y - center_y1) * (cloud->points[i].z - center_z1);
        zz1 += (cloud->points[i].z - center_z1) * (cloud->points[i].z - center_z1);

    }
    //大小为3*3的协方差矩阵
    Eigen::Matrix3d covMat1(3, 3);
    covMat1(0, 0) = xx1 / cld_sz_1;
    covMat1(0, 1) = covMat1(1, 0) = xy1 / cld_sz_1;
    covMat1(0, 2) = covMat1(2, 0) = xz1 / cld_sz_1;
    covMat1(1, 1) = yy1 / cld_sz_1;
    covMat1(1, 2) = covMat1(2, 1) = yz1 / cld_sz_1;
    covMat1(2, 2) = zz1 / cld_sz_1;

    //求特征值与特征向量
    Eigen::EigenSolver<Eigen::Matrix3d> es1(covMat1);
    Eigen::Matrix3d val1 = es1.pseudoEigenvalueMatrix();
    Eigen::Matrix3d vec1 = es1.pseudoEigenvectors();

    double* eigMatptr = val1.data();
    double* eigMatptrnew = new double[val1.size()];
    Eigen::Map<Eigen::Matrix3d>(eigMatptrnew, val1.rows(), val1.cols()) = val1;


    /*
    std::cout << val1.row(0) << "\n" << val1.row(1) << "\n" << val1.row(2) << std::endl;

    for (int i = 0; i < 9; ++i) {
        std::cout<< val1(i)  << std::endl;
    }

    */

    /*
    Eigen::RowVector3d EigenXYZ = vec1.row(0);
    Eigen::VectorXd a(1);
    a = EigenXYZ.col(0);
    */



    /*
    std::cout << "\n------" << dataname << "--------" << std::endl;
    std::cout << "     x  " << "     y  " << "     z  " << std::endl;
    std::cout << val1.row(0) << "\n" << val1.row(1) << "\n" << val1.row(2) << std::endl;
    std::cout << "\n------" << dataname << "--------" << std::endl;
    std::cout << "     x  " << "     y  " << "     z  " << std::endl;
    std::cout << vec1.row(0) << "\n" << vec1.row(1) << "\n" << vec1.row(2) << std::endl;
    */

    return val1;

    //return vec1;
}

int* sort(int num[], string dataname) {
    Eigen::Matrix3d matrix = eigen_value(dataname);

    double temp = 0;
    int tmpN = 0;

    std::vector<int> array;

    for (int i = 0; i < 9; i++) {
        for (int j = i; j < 9; j++) {
            if (matrix(j) > matrix(i)) {
                temp = matrix(j);
                tmpN = num[j];
                matrix(j) = matrix(i);
                num[j] = num[i];
                matrix(i) = temp;
                num[i] = tmpN;
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (num[i] == 4) {
            num[i] = 3;
        }
        if (num[i] == 8) {
            num[i] = 6;
        }
    }
    // 3的就是 Vector y (3'4'5) , 0 就是 vector x (0'1'2) , 6 就是 vector x (6'7'8) 
    return num;
}

struct node {
	int index;
	float curvature;
};

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
        newsearchPoint.y = (*cloud)[index].y + u * (vector(ptr[0]+1)) / 100;
        newsearchPoint.z = (*cloud)[index].z + u * (vector(ptr[0]+2)) / 100;

        int newK = 30;

        std::vector<int> newpointIdxNKNSearch(newK);
        std::vector<float> newpointNKNSquaredDistance(newK);

        float XPZa = vector(ptr[1]) * newsearchPoint.x;
        float XPZb = vector(ptr[1]+1) * newsearchPoint.y;
        float XPZc = vector(ptr[1]+2) * newsearchPoint.z;

        float XPYa = vector(ptr[2]) * newsearchPoint.x;
        float XPYb = vector(ptr[2]+1) * newsearchPoint.y;
        float XPYc = vector(ptr[2]+2) * newsearchPoint.z;

        if (newkdtree.nearestKSearch(newsearchPoint, newK, newpointIdxNKNSearch, newpointNKNSquaredDistance) > 0)
        {


            for (std::size_t i = 0; i < newpointIdxNKNSearch.size(); ++i) {

                float XOZa = vector(ptr[1]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float XOZb = vector(ptr[1]+1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float XOZc = vector(ptr[1]+2) * (*newcloud)[newpointIdxNKNSearch[i]].z;

                float XOYa = vector(ptr[2]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float XOYb = vector(ptr[2]+1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float XOYc = vector(ptr[2]+2) * (*newcloud)[newpointIdxNKNSearch[i]].z;

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
        for (int i = 0; i < 4; i++)
        {
            if (qindex.at(i) >= 0) {
                g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                    pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                    pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].z) / u1down;
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
        newsearchPoint.y = (*cloud)[index].y + u * (vector(ptr[1]+1)) / 100;
        newsearchPoint.z = (*cloud)[index].z + u * (vector(ptr[1]+2)) / 100;

        int newK = 30;

        std::vector<int> newpointIdxNKNSearch(newK);
        std::vector<float> newpointNKNSquaredDistance(newK);

        float YPZa = vector(ptr[0]) * newsearchPoint.x;
        float YPZb = vector(ptr[0]+1) * newsearchPoint.y;
        float YPZc = vector(ptr[0]+2) * newsearchPoint.z;

        float XPYa = vector(ptr[2]) * newsearchPoint.x;
        float XPYb = vector(ptr[2]+1) * newsearchPoint.y;
        float XPYc = vector(ptr[2]+2) * newsearchPoint.z;


        if (newkdtree.nearestKSearch(newsearchPoint, newK, newpointIdxNKNSearch, newpointNKNSquaredDistance) > 0)
        {

            for (std::size_t i = 0; i < newpointIdxNKNSearch.size(); ++i) {

                float YOZa = vector(ptr[0]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float YOZb = vector(ptr[0]+1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float YOZc = vector(ptr[0]+2) * (*newcloud)[newpointIdxNKNSearch[i]].z;

                float XOYa = vector(ptr[2]) * (*newcloud)[newpointIdxNKNSearch[i]].x;
                float XOYb = vector(ptr[2]+1) * (*newcloud)[newpointIdxNKNSearch[i]].y;
                float XOYc = vector(ptr[2]+2) * (*newcloud)[newpointIdxNKNSearch[i]].z;




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
        for (int i = 0; i < 4; i++)
        {
            if (qindex.at(i) >= 0) {
                g_ui += ((1 / (sqrt(pow(newsearchPoint.x - newcloud->points[qindex.at(i)].x, 2) +
                    pow(newsearchPoint.y - newcloud->points[qindex.at(i)].y, 2) +
                    pow(newsearchPoint.z - newcloud->points[qindex.at(i)].z, 2)))) * newcloud->points[qindex.at(i)].z) / u1down;
            }
        }

        g_kiy[u - 1] = g_ui;

    }
}

void reg(pcl::PointCloud<PointT>::Ptr target) {
    pcl::search::Search<PointT>::Ptr tree(new pcl::search::KdTree<PointT>);
    pcl::PointCloud <Normal>::Ptr normals(new pcl::PointCloud <Normal>);
    pcl::NormalEstimation<PointT, Normal> normal_estimator;
    normal_estimator.setSearchMethod(tree);
    normal_estimator.setInputCloud(target);
    normal_estimator.setKSearch(50);
    normal_estimator.compute(*normals);

    pcl::RegionGrowing<PointT, Normal> reg;
    reg.setMinClusterSize(100);
    reg.setMaxClusterSize(5000);
    reg.setSearchMethod(tree);
    reg.setNumberOfNeighbours(20);
    reg.setInputCloud(target);
    //reg.setIndices (indices);
    reg.setInputNormals(normals);
    reg.setSmoothnessThreshold(1.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold(0.1);

    std::vector <pcl::PointIndices> clusters;
    reg.extract(clusters);

    std::cout << "Number of clusters is equal to " << clusters.size() << std::endl;
    std::cout << "First cluster has " << clusters[0].indices.size() << " points." << std::endl;
    std::cout << "These are the indices of the points of the initial" <<
        std::endl << "cloud that belong to the first cluster:" << std::endl;
    int counter = 0;
    while (counter < clusters[0].indices.size())
    {
        std::cout << clusters[0].indices[counter] << ", ";
        counter++;
        if (counter % 10 == 0)
            std::cout << std::endl;
    }
    std::cout << std::endl;

    pcl::PointCloud <PointRGB>::Ptr colored_cloud = reg.getColoredCloud();

    pcl::io::savePCDFileASCII("reg_result1.pcd", *colored_cloud);
    std::cout << "saved!!" << std::endl;

    pcl::visualization::CloudViewer viewer("Cluster viewer");
    viewer.showCloud(colored_cloud);
    while (!viewer.wasStopped())
    {
    }
}

void printvector(vector<float> v) {
	for (int i = 0; i < v.size(); i++) {
		cout << v[i] << endl;
	}
}

void Merge(std::vector<node>& A, int i, int m, int j) { //A用pointer表示
	std::vector<node> Leftsub(A.begin() + i, A.begin() + m + 1);
	std::vector<node> Rightsub(A.begin() + m + 1, A.begin() + j + 1);
	node max;
	max.index = -1; max.curvature = Min;
	Leftsub.push_back(max); //兩個sublist後面加一個無限大
	Rightsub.push_back(max); //兩個sublist後面加一個無限大

	int L = 0;
	int R = 0;
	for (int f = i; f <= j; f++) {
		if (L>=Leftsub.size()-1 && R>=Rightsub.size()-1) break;
		else {
            if (R >= Rightsub.size() - 1 || L >= Leftsub.size() - 1) {
                if (R >= Rightsub.size() - 1) {
                    A[f] = Leftsub[L];
                    L++;
                }
                else {
                    A[f] = Rightsub[R];
                    R++;
                }
            }
            else if((fabsf(Leftsub[L].curvature) > fabsf(Rightsub[R].curvature)) || (fabsf(Rightsub[R].curvature) > fabsf(Leftsub[L].curvature))) { //由大到小
                if (fabsf(Leftsub[L].curvature) > fabsf(Rightsub[R].curvature)) {
                    A[f] = Leftsub[L];
                    L++;
                }
                else {
                    A[f] = Rightsub[R];
                    R++;
                }
			}
			else {
                if (L >= Leftsub.size()-1) {
                    A[f] = Rightsub[R];
                    R++;
                }
                else {
                    A[f] = Leftsub[L];
				    L++;
                }
			}
		}
	}
}

void mergesort(std::vector<node>& A, int i, int j) {
	if (i < j)
	{
		int m = (i + j) / 2;
		mergesort(A, i, m); //devide &conquer
		mergesort(A, m + 1, j); //devide &conquer
		Merge(A, i, m, j); //combine
	}
}

void PrintArray(std::vector<node>& array) {
	for (int i = 0; i < array.size(); i++) {
		std::cout << array[i].index << " : " << array[i].curvature << std::endl;
	}
	std::cout << std::endl;
}

void transform_Matrix(string basedata, string transformresult, int Xangle, int Yangle, int Zangle) {
    typedef pcl::PointXYZ PointT;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::io::loadPCDFile(basedata + ".pcd", *cloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
    //pcl::io::loadPCDFile("A0022_tc_test8.pcd", *source_cloud);
    Eigen::Matrix4f transform_1 = Eigen::Matrix4f::Identity();

    float sa = sin(Xangle * PI / 180);
    float ca = cos(Xangle * PI / 180);

    float sb = sin(Yangle * PI / 180);
    float cb = cos(Yangle * PI / 180);

    float sc = sin(Zangle * PI / 180);
    float cc = cos(Zangle * PI / 180);

    // x y z 旋轉
    transform_1(0, 0) = cb * cc;                       transform_1(0, 1) = cb * sc;                     transform_1(0, 2) = -sb;          transform_1(0, 3) = 0;

    transform_1(1, 0) = -ca * sc + sa * sb * cc;       transform_1(1, 1) = ca * cc + sa * sb * sc;      transform_1(1, 2) = sa * cb;      transform_1(1, 3) = 0;

    transform_1(2, 0) = sa * sc + ca * sb * cc;        transform_1(2, 1) = -sa * cc + ca * sb * sc;     transform_1(2, 2) = ca * cb;      transform_1(2, 3) = 0;

    transform_1(3, 0) = 0.0;                           transform_1(3, 1) = 0.0;                         transform_1(3, 2) = 0.0;          transform_1(3, 3) = 1.0;

    pcl::transformPointCloud(*cloud, *transformed_cloud, transform_1);
    pcl::io::savePCDFile<pcl::PointXYZ>(transformresult + ".pcd", *transformed_cloud);
}

Eigen::Matrix3d eigenv_v(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) {
    int cld_sz_1 = cloud->size();
    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    normals->resize(cld_sz_1);

    //计算中心点坐标
    double center_x1 = 0, center_y1 = 0, center_z1 = 0;
    for (int i = 0; i < cld_sz_1; i++) {
        center_x1 += cloud->points[i].x;
        center_y1 += cloud->points[i].y;
        center_z1 += cloud->points[i].z;

    }
    center_x1 /= cld_sz_1;
    center_y1 /= cld_sz_1;
    center_z1 /= cld_sz_1;
    //计算协方差矩阵
    double xx1 = 0, xy1 = 0, xz1 = 0, yy1 = 0, yz1 = 0, zz1 = 0;
    for (int i = 0; i < cld_sz_1; i++) {
        xx1 += (cloud->points[i].x - center_x1) * (cloud->points[i].x - center_x1);
        xy1 += (cloud->points[i].x - center_x1) * (cloud->points[i].y - center_y1);
        xz1 += (cloud->points[i].x - center_x1) * (cloud->points[i].z - center_z1);
        yy1 += (cloud->points[i].y - center_y1) * (cloud->points[i].y - center_y1);
        yz1 += (cloud->points[i].y - center_y1) * (cloud->points[i].z - center_z1);
        zz1 += (cloud->points[i].z - center_z1) * (cloud->points[i].z - center_z1);

    }
    //大小为3*3的协方差矩阵
    Eigen::Matrix3d covMat1(3, 3);
    covMat1(0, 0) = xx1 / cld_sz_1;
    covMat1(0, 1) = covMat1(1, 0) = xy1 / cld_sz_1;
    covMat1(0, 2) = covMat1(2, 0) = xz1 / cld_sz_1;
    covMat1(1, 1) = yy1 / cld_sz_1;
    covMat1(1, 2) = covMat1(2, 1) = yz1 / cld_sz_1;
    covMat1(2, 2) = zz1 / cld_sz_1;

    //求特征值与特征向量
    Eigen::EigenSolver<Eigen::Matrix3d> es1(covMat1);
    Eigen::Matrix3d val1 = es1.pseudoEigenvalueMatrix();
    Eigen::Matrix3d vec1 = es1.pseudoEigenvectors();

    double* eigMatptr = vec1.data();
    double* eigMatptrnew = new double[vec1.size()];
    Eigen::Map<Eigen::Matrix3d>(eigMatptrnew, vec1.rows(), vec1.cols()) = vec1;

    for (int i = 0; i < 9; ++i) {
        //std::cout<< vec1(i) + 1 << std::endl;
    }

    /*
    Eigen::RowVector3d EigenXYZ = vec1.row(0);
    Eigen::VectorXd a(1);
    a = EigenXYZ.col(0);
    */

    /*
    std::cout << "\n------" << dataname << "--------" << std::endl;
    std::cout << "     x  " << "     y  " << "     z  " << std::endl;
    std::cout << val1.row(0) << "\n" << val1.row(1) << "\n" << val1.row(2) << std::endl;
    std::cout << "\n------" << dataname << "--------" << std::endl;
    std::cout << "     x  " << "     y  " << "     z  " << std::endl;
    std::cout << vec1.row(0) << "\n" << vec1.row(1) << "\n" << vec1.row(2) << std::endl;
    */

    return vec1;
}

void simp(/*pcl::PointCloud<PointT>::Ptr input*/)  //主程式
{
    std::vector<node> h_mean; //主要儲存vector
    std::string filename = "S0001A0022_XYZ";
    //load the file
    pcl::PointCloud<PointT>::Ptr input(new pcl::PointCloud<PointT>);


    if (pcl::io::loadPCDFile(filename+".pcd", *input) == -1) {
        PCL_ERROR("Can not read the file.");
    }

    int num[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    //std::cout << eigen_value(filename) << std::endl;

    int* ptr = sort(num, filename);

    for (int i = 0; i < 3; i++) {
        std::cout << ptr[i] << std::endl;
    }

    Eigen::Matrix3d vector_eigen = eigenv_v(input);

    //build the kdtree
    pcl::KdTreeFLANN<PointT> kdtree;
    kdtree.setInputCloud(input);

    cout << "start sampling" << endl;
    for (int i = 0; i < input->size(); i++) {

        auto start = std::chrono::system_clock::now();
        // Some computation here


        PointT searchPoint;
        searchPoint = input->points[i];
        pcl::PointCloud<PointT>::Ptr newcloud(new pcl::PointCloud<PointT>);
        int K = 31;
        std::vector<int> pointIdxNKNSearch(K);
        std::vector<float> pointNKNSquaredDistance(K);


        if (kdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
        {


            for (std::size_t i = 1; i < pointIdxNKNSearch.size(); ++i) {


                newcloud->push_back((*input)[pointIdxNKNSearch[i]]);

            }

        }
        pcl::KdTreeFLANN<PointT> newkdtree;


        newkdtree.setInputCloud(newcloud);

        //auto end = std::chrono::system_clock::now();

        //std::chrono::duration<double> elapsed_seconds = end - start;

        //std::cout << elapsed_seconds.count() << std::endl;


        kix(i, input, newcloud, newkdtree, vector_eigen, ptr);
        kiy(i, input, newcloud, newkdtree, vector_eigen, ptr);


        /*
    *fftw_complex 是FFTW自定义的复数类
    *引入<complex>则会使用STL的复数类
    */

        int N = 8;
        double Ki_first = 0;
        double Ki_second = 0;
        double H_i;

        node point;


        Ki_first = k_i(g_kix, N);

        Ki_second = k_i(g_kiy, N);

        H_i = (Ki_first + Ki_second) / 2;


        point.index = i;
        point.curvature = H_i;

        h_mean.push_back(point);

        //cout << "h_mean index  " << h_mean.at(i).index<< "h_mean curvature  " << h_mean.at(i).curvature << endl;

    }
    cout << "finish the curvature" << endl;

    float Mean = 0;
    cout << "start mean" << endl;
    for (int c = 0; c < h_mean.size(); c++)
        Mean += fabsf(h_mean.at(c).curvature);
    Mean = Mean / input->size();
    cout << "finish mean" << endl;

    vector<int> indices;
    //vector<int> record;
    vector<float> dist;
    //vector<float> curv;
    vector<int> checked;
    vector<int> storage_index;

    cout << "start merging" << endl;
    mergesort(h_mean, 0, h_mean.size() - 1);
    cout << "finish merging" << endl;

    //cout << "test stop" << endl;

    for (int i = 0; i < input->size(); i++)
        checked.push_back(0);

    pcl::search::KdTree<PointT>::Ptr Kdtree(new pcl::search::KdTree<PointT>);
    Kdtree->setInputCloud(input);

    //int k = 0;

    float alpha = 0.14;
    cout << "start marking" << endl;
    pcl::PointXYZ searchPoint;
    for (int k = 0; k < h_mean.size(); k++) {
        if (checked[h_mean[k].index] == 0) {
            searchPoint.x = input->points[h_mean[k].index].x;
            searchPoint.y = input->points[h_mean[k].index].y;
            searchPoint.z = input->points[h_mean[k].index].z;

            //cout << searchPoint.x << " " << searchPoint.y << " " << searchPoint.z << endl;

            float radius = alpha * Mean / fabsf(h_mean[k].curvature);
            Kdtree->radiusSearch(searchPoint, radius, indices, dist);

            //cout << "cloud size after search:" << input->size() << endl;

            for (int j = 0; j < indices.size(); j++) {
                if (checked[indices[j]] != 2)
                    checked[indices[j]] = 1;
            }
            checked[h_mean[k].index] = 2;
            indices.clear();
        }
    }
    cout << "finish marking" << endl;


    for (int tap = 0; tap < checked.size(); tap++) {
        if (checked[tap] == 2)
            storage_index.push_back(tap);
    }

    //cout << "stop" << endl;
    int counted = storage_index.size();
    //cout << "remain points:" << counted << endl;
    //cout << storage_index.size() << endl;

    pcl::PointCloud<pcl::PointXYZ> result;
    result.width = counted;
    result.height = 1;
    result.points.resize(result.width * result.height);
    int trace = 0;
    for (auto& point : result)
    {
        point.x = input->points[storage_index[trace]].x;
        point.y = input->points[storage_index[trace]].y;
        point.z = input->points[storage_index[trace]].z;
        trace++;
    }
    cout << "ready to save" << endl;
    pcl::io::savePCDFileASCII(filename + "simp001.pcd", result);
}

int main() {
    simp();
    
}