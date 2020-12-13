/*
 * This was Hough transform algorithm but works bad
 */
#include <iostream>
//#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstring>
#include <pangolin/pangolin.h>
#include <ctime>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/passthrough.h>

using namespace std;
#define PI 3.141592653

//class Point
//{
//public:
//    double x, y, z;
//    int planID;
//    Point(double ix,double iy,double iz) :
//            x(ix), y(iy), z(iz){planID=-1;}
//
//    Point operator-(const Point& pt) const
//    {
//        return Point(x - pt.x, y - pt.y, z - pt.z);
//    }
//};
//
//class Plane {
//public:
//    double A, B, C, D;
//    double phi, theta, dis;
//    int count;
//    vector<int> pointList;
//    int mergeID;
//    Plane(double Ain, double Bin, double Cin, double Din) : A(Ain), B(Bin), C(Cin), D(Din) {mergeID=-1;count=-1;};
//};

int main() {
    std::cout << "Hello, World!" << std::endl;

//    ///define pcl kdtree stuff
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
//    // Generate pointcloud data
//    cloud->width = 100000;
//    cloud->height = 1;
//    for(size_t i=0; i<cloud->points.size();i++)
//    {
//        cloud->points[i].x = 1024.0f * rand()/(RAND_MAX+1.0f);
//        cloud->points[i].y = 1024.0f * rand()/(RAND_MAX+1.0f);
//        cloud->points[i].z = 1024.0f * rand()/(RAND_MAX+1.0f);
//    }
    cloud->points.resize (200000);

///load KITII data
    //Load data from bin
    string dataFile = "../0000000000.bin";
    int32_t num = 5000000;
    float *data = (float *) malloc(num * sizeof(float));
    float *px = data + 0;
    float *py = data + 1;
    float *pz = data + 2;
    float *pr = data + 3;
    FILE *fstream;
    fstream = fopen(dataFile.c_str(), "rb");
    num = fread(data, sizeof(float), num, fstream) / 4;
    int actualNum = 0;
    for (int i = 0; i < num; i++) {
        if (*px > 0 && *px < 25 && abs(*py) < 6) {
            cloud->points[actualNum].x = *px;
            cloud->points[actualNum].y = *py;
            cloud->points[actualNum].z = *pz;
            actualNum++;
        }
        px += 4;
        py += 4;
        pz += 4;
        pr += 4;
    }
    fclose(fstream);
    cout << "actual loaded " << actualNum << " points" << endl;
    cloud->points.resize (actualNum);


    //create instance of kdtreeflann
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud);
    //set search point
    pcl::PointXYZ searchPoint;
    searchPoint.x = 1024.0f * rand () / (RAND_MAX + 1.0f);
    searchPoint.y = 1024.0f * rand () / (RAND_MAX + 1.0f);
    searchPoint.z = 1024.0f * rand () / (RAND_MAX + 1.0f);

    //K nearest neighbour search
    int K =10;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    std::cout << "K nearest neighbor search at (" << searchPoint.x
              << " " << searchPoint.y
              << " " << searchPoint.z
              << ") with K=" << K << std::endl;

    if(kdtree.nearestKSearch(searchPoint,K, pointIdxNKNSearch, pointNKNSquaredDistance)>0)
    {
        for (size_t i = 0; i < pointIdxNKNSearch.size (); ++i)
            std::cout << "    "  <<   cloud->points[ pointIdxNKNSearch[i] ].x
                      << " " << cloud->points[ pointIdxNKNSearch[i] ].y
                      << " " << cloud->points[ pointIdxNKNSearch[i] ].z
                      << " (squared distance: " << pointNKNSquaredDistance[i] << ")" << std::endl;
    }

    // 半径 R内近邻搜索方法
    std::vector<int> pointIdxRadiusSearch;           //存储近邻索引
    std::vector<float> pointRadiusSquaredDistance;   //存储近邻对应距离的平方

    float radius = 256.0f * rand () / (RAND_MAX + 1.0f);   //随机的生成某一半径
    //打印输出
    std::cout << "Neighbors within radius search at (" << searchPoint.x
              << " " << searchPoint.y
              << " " << searchPoint.z
              << ") with radius=" << radius << std::endl;


    if ( kdtree.radiusSearch (searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0 )  //执行半径R内近邻搜索方法
    {
        for (size_t i = 0; i < pointIdxRadiusSearch.size (); ++i)
            std::cout << "    "  <<   cloud->points[ pointIdxRadiusSearch[i] ].x
                      << " " << cloud->points[ pointIdxRadiusSearch[i] ].y
                      << " " << cloud->points[ pointIdxRadiusSearch[i] ].z
                      << " (squared distance: " << pointRadiusSquaredDistance[i] << ")" << std::endl;
    }

    ///segment
    //estimating normals for each point
    pcl::search::Search<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
    pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator;
    normal_estimator.setSearchMethod (tree);
    normal_estimator.setInputCloud (cloud);
    normal_estimator.setKSearch (50);
    normal_estimator.compute (*normals);

//    pcl::IndicesPtr indices (new std::vector <int>);
//    pcl::PassThrough<pcl::PointXYZ> pass;
//    pass.setInputCloud (cloud);
//    pass.setFilterFieldName ("z");
//    pass.setFilterLimits (0.0, 1.0);
//    pass.filter (*indices);

    pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
    reg.setMinClusterSize (1000);
    reg.setMaxClusterSize (1000000);
    reg.setSearchMethod (tree);
    reg.setNumberOfNeighbours (30);
    reg.setInputCloud (cloud);
    //reg.setIndices (indices);
    reg.setInputNormals (normals);
    reg.setSmoothnessThreshold (3.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold (1.0);

    std::vector <pcl::PointIndices> clusters;
    reg.extract (clusters);

//    std::cout << "Number of clusters is equal to " << clusters.size () << std::endl;
//    std::cout << "First cluster has " << clusters[0].indices.size () << " points." << std::endl;
//    std::cout << "These are the indices of the points of the initial" <<
//              std::endl << "cloud that belong to the first cluster:" << std::endl;
//    int counter = 0;
//    while (counter < clusters[0].indices.size ())
//    {
//        //std::cout << clusters[0].indices[counter] << ", ";
//        counter++;
//        if (counter % 500 == 0)
//        {
//            std::cout << clusters[0].indices[counter] << ", ";
//            //std::cout<<cloud->points[counter].x<<" "<<cloud->points[counter].y<<" "<<cloud->points[counter].z<<" ";
//            std::cout << std::endl;
//        }
//    }
//    std::cout << std::endl;

    int showCount = 0;
    for(int ci = 0; ci<clusters.size();ci++)
    {
        cout<<"cluster "<<ci<<" has "<<clusters[ci].indices.size()<<" points :"<<endl;
        int pi = 0;
        while(pi<clusters[ci].indices.size())
        {
            int Pindex = clusters[ci].indices[pi];
            if(pi % 500 == 0 && showCount < 10)
            {
                cout<<cloud->points[Pindex].x<<" "<<cloud->points[Pindex].y<<" "<<cloud->points[Pindex].z<<endl;
            }
            pi++;
        }
    }


    pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud ();
    pcl::visualization::CloudViewer viewer ("Cluster viewer");
    viewer.showCloud(colored_cloud);
    while (!viewer.wasStopped ())
    {
    }

    return (0);


//
//
//
//
//    ///show
//    pangolin::CreateWindowAndBind("Main", 640, 480);
//    glEnable(GL_DEPTH_TEST);
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//
//    pangolin::OpenGlRenderState s_cam(
//            pangolin::ProjectionMatrix(640, 480, 420, 420, 320, 320, 0.2, 100),
//            pangolin::ModelViewLookAt(2, 0, 2, 0, 0, 0, pangolin::AxisY)
//    );
//
//    pangolin::Handler3D handler(s_cam);
//    pangolin::View &d_cam = pangolin::CreateDisplay()
//            .SetBounds(0.0, 1.0, 0.0, 1.0, -640.0f / 480.0f)
//            .SetHandler(&handler);
//
//    while (!pangolin::ShouldQuit()) {
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//        d_cam.Activate(s_cam);
//
//        //需要绘制的东西写在这里
//        for (int pi = 0; pi < laserPoints.size(); pi++) {
//            glBegin(GL_POINTS);
//            switch (laserPoints[pi].planID){
//                case -1 :
//                    glColor3f(0.8, 0.8, 0.8);
//                    break;
//                case 0 :
//                    glColor3f(1, 0, 0);
//                    break;
//                case 1 :
//                    glColor3f(0, 1, 0);
//                    break;
//                case 2 :
//                    glColor3f(0, 0, 1);
//                    break;
//                case 3 :
//                    glColor3f(1, 1, 0);
//                    break;
//                case 4 :
//                    glColor3f(0, 1, 1);
//                    break;
//            }
//            glVertex3d(laserPoints[pi].x, laserPoints[pi].y, laserPoints[pi].z);
//            glEnd();//点设置的结束
//        }
//        //box of max plane
//        glBegin(GL_LINES);
//        glLineWidth(2.0);
//        glColor3f(1,1,1);
//        glVertex3d(selPoint.x,selPoint.y,selPoint.z);
//        glVertex3d(selPoint.x+maxA,selPoint.y+maxB,selPoint.z+maxC);
//        glEnd();
//        //coordinate frame
//        glBegin(GL_LINES);
//        glLineWidth(5.0);
//        glColor3f(1,0,1);
//        glVertex3d(0.1,0,0);
//        glVertex3d(0,0,0);
//        glColor3f(0,1,0);
//        glVertex3d(0,0.1,0);
//        glVertex3d(0,0,0);
//        glColor3f(0,0,1);
//        glVertex3d(0,0,0.1);
//        glVertex3d(0,0,0);
//        glEnd();
//        //box of unit plane
//        glBegin(GL_LINES);
//        glLineWidth(10.0);
//        glColor3f(1,1,1);
//        glVertex3d(uniplanecorners[0].x,uniplanecorners[0].y,uniplanecorners[0].z);
//        glVertex3d(uniplanecorners[1].x,uniplanecorners[1].y,uniplanecorners[1].z);
//        glVertex3d(uniplanecorners[0].x,uniplanecorners[0].y,uniplanecorners[0].z);
//        glVertex3d(uniplanecorners[2].x,uniplanecorners[2].y,uniplanecorners[2].z);
//        glVertex3d(uniplanecorners[1].x,uniplanecorners[1].y,uniplanecorners[1].z);
//        glVertex3d(uniplanecorners[3].x,uniplanecorners[3].y,uniplanecorners[3].z);
//        glVertex3d(uniplanecorners[2].x,uniplanecorners[2].y,uniplanecorners[2].z);
//        glVertex3d(uniplanecorners[3].x,uniplanecorners[3].y,uniplanecorners[3].z);
//        glColor3f(1,1,1);
//        glVertex3d(0,0,0);
//        glVertex3d(0,0,0.05);
//        glEnd();
//        //
//
//        pangolin::FinishFrame();
//    }

    return 0;
}
