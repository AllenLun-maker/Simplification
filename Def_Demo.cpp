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
#include <pcl/visualization/cloud_viewer.h>


//#include "Measure.h"

using namespace std;

typedef pcl::PointXYZ PointT;
typedef pcl::PointXYZRGB PointRGB;
typedef pcl::Normal Normal;


string inputcloud;

class Show {
public :
    void Visualize() {
            pcl::PointCloud<pcl::PointXYZ>::Ptr input(new pcl::PointCloud<pcl::PointXYZ>);
            if (pcl::io::loadPCDFile(inputcloud + ".pcd", *input) == -1) {
                PCL_ERROR("Can not read the file.");
            }
            pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("VIEWER"));
            viewer->setBackgroundColor(0, 0, 0);
            viewer->addPointCloud<pcl::PointXYZ>(input, "cloud");
            while (!viewer->wasStopped())
            {
                viewer->spinOnce(1000);
                //boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            }
            std::cout << "Visualize complete";
        }
};




