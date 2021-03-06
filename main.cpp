/*
 * This was Region growing first (0.3 sec)
 * then RANSAC to fit plane for each cluster
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
#include <pcl/segmentation/region_growing.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/ModelCoefficients.h>

using namespace std;
double PI = 3.1415926;

class Point {
public:
    double x, y, z;
    int planID;

    Point(double ix, double iy, double iz) :
            x(ix), y(iy), z(iz) { planID = -1; }

    Point operator-(const Point &pt) const {
        return Point(x - pt.x, y - pt.y, z - pt.z);
    }
};

class Plane {
public:
    double A, B, C, D;
    double phi, theta, dis;
    int count;
    vector<int> pointList;
    int mergeID;

    Plane(double Ain, double Bin, double Cin, double Din, Point centreP)
            : A(Ain), B(Bin), C(Cin), D(Din), centreP(centreP) {
        mergeID = -1;
        count = -1;
    };

    Plane(Point centrePin) : centreP(centrePin) {};
    Point centreP;
};

void RANSAC_plan(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, Plane &foundPlane)
{
    //pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = inputCloud.makeShared();
    pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
    //create the segmentation object
    pcl::SACSegmentation<pcl::PointXYZ> seg;
    //Optional
    seg.setOptimizeCoefficients(true);
    //Mandatory
    seg.setModelType(pcl::SACMODEL_PLANE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setDistanceThreshold(0.01);

    seg.setInputCloud(cloud);
    seg.segment(*inliers, *coefficients);
    if(inliers->indices.size()==0)
    {
        PCL_ERROR("inlier number 0");
        return;
    }

    std::cout<<"Model coefficients: "<<coefficients->values[0]<<" "<<coefficients->values[1]<<" "<<coefficients->values[2]<<" "<<coefficients->values[3]<<endl;
    std::cout<<"Model inlier number: "<<inliers->indices.size()<<endl;
    foundPlane.A = coefficients->values[0];foundPlane.B = coefficients->values[1];
    foundPlane.C = coefficients->values[2];foundPlane.D = coefficients->values[3];
    //TODO add point index to plane
    double sumX =0,sumY=0,sumZ=0;
    for(int i = 0; i < inliers->indices.size();i++)
    {
        sumX += cloud->points[inliers->indices[i]].x;
        sumY += cloud->points[inliers->indices[i]].y;
        sumZ += cloud->points[inliers->indices[i]].z;
    }
    foundPlane.centreP = Point(sumX/inliers->indices.size(), sumY/inliers->indices.size(),sumZ/inliers->indices.size());
}

int main() {

    ///Load data
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    cloud->points.resize(1000000);
    std::string folderName = "/home/xin/Downloads/DATASET/KITTI/RAW/2011_09_26_drive_0046/2011_09_26_drive_0046_sync/2011_09_26/2011_09_26_drive_0046_sync/velodyne_points/data/";
    std::string fronter1 = "000000000";
    std::string fronter2 = "00000000";
    std::string fronter3 = "0000000";
    std::string ender = ".bin";
    std::string fileAddress = "";
    for (int i = 0; i < 125; i++) {
        stringstream ss;
        ss << i;
        if (i < 10)
            fileAddress = folderName + fronter1 + ss.str() + ender;
        else if (i < 100)
            fileAddress = folderName + fronter2 + ss.str() + ender;
        else if (i >= 100)
            fileAddress = folderName + fronter3 + ss.str() + ender;
        cout << fileAddress << endl;

        ///load KITTI data
        int32_t num = 1000000;
        float *data = (float *) malloc(num * sizeof(float));
        float *px = data + 0;
        float *py = data + 1;
        float *pz = data + 2;
        float *pr = data + 3;
        FILE *fstream;
        fstream = fopen(fileAddress.c_str(), "rb");
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
        cloud->points.resize(actualNum);

        ///estimating normals for each point
        pcl::search::Search<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
        pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
        pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator;
        normal_estimator.setSearchMethod(tree);
        normal_estimator.setInputCloud(cloud);
        normal_estimator.setKSearch(50);
        normal_estimator.compute(*normals);

        //    pcl::IndicesPtr indices (new std::vector <int>);
        //    pcl::PassThrough<pcl::PointXYZ> pass;
        //    pass.setInputCloud (cloud);
        //    pass.setFilterFieldName ("z");
        //    pass.setFilterLimits (0.0, 1.0);
        //    pass.filter (*indices);

        ///region growing
        pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
        reg.setMinClusterSize(1000);
        reg.setMaxClusterSize(1000000);
        reg.setSearchMethod(tree);
        reg.setNumberOfNeighbours(100); //too little will wrong time error.
        reg.setInputCloud(cloud);
        //reg.setIndices (indices);
        reg.setInputNormals(normals);
        reg.setSmoothnessThreshold(3.0 / 180.0 * M_PI);
        reg.setCurvatureThreshold(2.0);

        clock_t startTime = clock();
        std::vector<pcl::PointIndices> clusters;
        reg.extract(clusters);
        clock_t endTime = clock();
        double timeUsed = double(endTime - startTime) / CLOCKS_PER_SEC;
        cout << "Region Growing time used " << timeUsed << " sec" << endl;

        ///PCL show
//        pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
//        pcl::visualization::CloudViewer viewer("Cluster viewer");
//        viewer.showCloud(colored_cloud);
//        while (!viewer.wasStopped()) {
//        }

        ///plane finding
        vector<Plane> planes;
        for (int ci = 0; ci < clusters.size(); ci++) {
            pcl::PointCloud<pcl::PointXYZ>::Ptr thisCloud(new pcl::PointCloud<pcl::PointXYZ>);
            thisCloud->points.resize(clusters[ci].indices.size());
            thisCloud->height=1;thisCloud->width=clusters[ci].indices.size();
            for (int index = 0; index < clusters[ci].indices.size(); index++) {
                thisCloud->points[index].x = cloud->points[clusters[ci].indices[index]].x;
                thisCloud->points[index].y = cloud->points[clusters[ci].indices[index]].y;
                thisCloud->points[index].z = cloud->points[clusters[ci].indices[index]].z;
            }
            Plane foundPlane(Point(0,0,0));
            clock_t startTime_ = clock();
            RANSAC_plan(thisCloud, foundPlane);
            clock_t endTime_ = clock();
            double timeUsed_ = double(endTime_-startTime_)/CLOCKS_PER_SEC;
            cout<<"RANSAC time used "<<timeUsed_<<endl;
            planes.push_back(foundPlane);
        }

        ///show
        pangolin::CreateWindowAndBind("Main", 640, 480);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        pangolin::OpenGlRenderState s_cam(
                pangolin::ProjectionMatrix(640, 480, 420, 420, 320, 320, 0.2, 100),
                pangolin::ModelViewLookAt(2, 0, 2, 0, 0, 0, pangolin::AxisY)
        );

        pangolin::Handler3D handler(s_cam);
        pangolin::View &d_cam = pangolin::CreateDisplay()
                .SetBounds(0.0, 1.0, 0.0, 1.0, -640.0f / 480.0f)
                .SetHandler(&handler);
        while (!pangolin::ShouldQuit()) {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            d_cam.Activate(s_cam);

            //?????????????????????????????????
            //coordinate frame
            glBegin(GL_LINES);
            glLineWidth(5.0);
            glColor3f(1, 0, 1);
            glVertex3d(0.1, 0, 0);
            glVertex3d(0, 0, 0);
            glColor3f(0, 1, 0);
            glVertex3d(0, 0.1, 0);
            glVertex3d(0, 0, 0);
            glColor3f(0, 0, 1);
            glVertex3d(0, 0, 0.1);
            glVertex3d(0, 0, 0);
            glEnd();
            //points
            for (int ci = 0; ci < clusters.size(); ci++) {
                for (int pi = 0; pi < clusters[ci].indices.size(); pi++) {
                    glBegin(GL_POINTS);
                    switch (ci) {
                        case 0 :
                            glColor3f(1, 0, 0);
                            break;
                        case 1 :
                            glColor3f(0, 1, 0);
                            break;
                        case 2 :
                            glColor3f(0, 0, 1);
                            break;
                        case 3 :
                            glColor3f(1, 1, 0);
                            break;
                        case 4 :
                            glColor3f(1, 0, 1);
                            break;
                        case 5 :
                            glColor3f(0, 1, 1);
                            break;
                    }
                    glVertex3d(cloud->points[clusters[ci].indices[pi]].x, cloud->points[clusters[ci].indices[pi]].y, cloud->points[clusters[ci].indices[pi]].z);
                    glEnd();//??????????????????
                }
            }

            //normals of planes
            glBegin(GL_LINES);
            glLineWidth(2.0);
            glColor3f(1, 1, 1);
            for(int pi = 0; pi<planes.size();pi++)
            {
                Point cp = planes[pi].centreP;
                Point np = Point(planes[pi].A,planes[pi].B,planes[pi].C);
                np.x = np.x + cp.x; np.y = np.y + cp.y; np.z = np.z + cp.z;
                glVertex3d(cp.x,cp.y,cp.z);
                glVertex3d(np.x,np.y,np.z);
            }

            glEnd();

            pangolin::FinishFrame();
        }


    }
    return 0;
}
