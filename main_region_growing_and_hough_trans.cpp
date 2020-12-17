/*
 * This was Region growing first (0.3 sec)
 * then Hough trans for finding plane in each cluster
 * time consumes 0.2(each small cluster) 1.8 (larger cluster)
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
#include <pcl/filters/passthrough.h>

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

class boundingbox {
public:
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;

public:
    double diag() const {
        double dx = x_max - x_min;
        double dy = y_max - y_min;
        double dz = z_max - z_min;

        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    boundingbox() :
            x_min(std::numeric_limits<double>::max()),
            y_min(std::numeric_limits<double>::max()),
            z_min(std::numeric_limits<double>::max()),
            x_max(-std::numeric_limits<double>::max()),
            y_max(-std::numeric_limits<double>::max()),
            z_max(-std::numeric_limits<double>::max()) {}

    Point center() const {
        return Point((x_max + x_min) / 2.0, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0);
    }
};

void calcboundbox(const std::vector<Point> &input, boundingbox &box) {
    for (int i = 0, n = input.size(); i < n; ++i) {
        auto point = input[i];
        if (point.x < box.x_min)
            box.x_min = point.x;
        if (point.y < box.y_min)
            box.y_min = point.y;
        if (point.z < box.z_min)
            box.z_min = point.z;
        if (point.x > box.x_max)
            box.x_max = point.x;
        if (point.y > box.y_max)
            box.y_max = point.y;
        if (point.z > box.z_max)
            box.z_max = point.z;
    }
}

void HoughTransform(const vector<Point> &input, double &A, double &B, double &C, double &D, Plane &foundPlane, int mode) {
    int n = input.size();
    if (n < 3)
        return;

    double theta_start = 0, theta_end = PI;
    double phi_start = 0, phi_end = PI;
    //double phi_start = 0.25*PI, phi_end = 0.75*PI;
    double anglestep = PI / 90, disstep = 0.05;
    //floor
    if(mode ==1)
    {
        theta_start = PI / 2 - 0.2; theta_end = PI / 2 + 0.2;
    }
    if(mode ==2)
    {
        phi_start = PI/2-0.2, phi_end = PI/2+0.2;
    }

    boundingbox box;
    calcboundbox(input, box);
    double d_start = -box.diag() / 2.0, d_end = box.diag() / 2.0;

    int thetas = ceil((theta_end - theta_start) / anglestep);
    int phis = ceil((phi_end - phi_start) / anglestep);
    int dises = ceil(box.diag() / disstep);

    int ***cube = new int **[thetas];
    for (int i = 0; i < thetas; ++i) {
        cube[i] = new int *[phis];
        for (int j = 0; j < phis; ++j) {
            cube[i][j] = new int[dises];
            memset(cube[i][j], 0, sizeof(int) * dises);
        }
    }

    //cos(theta)sin(phi)X+sin(theta)sin(phi)Y+cos(phi)Z = D
    Point ptCenter = box.center();
    for (int i = 0; i < n; ++i) {
        const Point &ptOrigin = input[i];
        Point point = ptOrigin - ptCenter;

        double theta = theta_start;
        for (int j = 0; j < thetas; ++j) {
            int **row = cube[j];
            double phi = phi_start;
            for (int k = 0; k < phis; ++k) {
                int *col = row[k];

                double sinphi = sin(phi);
                double d = cos(theta) * sinphi * point.x + sin(theta) * sinphi * point.y + cos(phi) * point.z;
                int d_index = floor((d - d_start) / disstep);

//                if (i % 1000 == 0)
//                    cout << "dealing with point " << i << " " << " theta " << theta << " index " << j << " phi " << phi
//                         << " index " << k << " dis " << d << " index " << d_index << endl;
                ++(col[d_index]);

                phi += anglestep;

                if (phi > phi_end)
                    break;
            }
            theta += anglestep;
            if (theta > theta_end)
                break;
        }
    }//all points

    int buf = 1;
    int maxcount = 0;
    int xmax, ymax, zmax;
    int countThres = 80000;
    for (int i = 0; i < thetas; ++i)
        for (int j = 0; j < phis; ++j)
            for (int k = buf; k < dises - buf; ++k) {
                int count = 0;
                for (int x = i - buf; x <= i + buf; ++x)
                    for (int y = j - buf; y <= j + buf; ++y)
                        for (int z = k - buf; z <= k + buf; ++z) {
                            count += cube[x < 0 ? x + thetas : x % thetas][y < 0 ? y + phis : y % phis][z];
                        }
                if (count > maxcount) {
                    xmax = i;
                    ymax = j;
                    zmax = k;
                    maxcount = count;
                    //cout<<"new max count "<<maxcount<<endl;
                }
            }

    double theta = theta_start + xmax * anglestep;
    double phi = phi_start + ymax * anglestep;
    double d = d_start + zmax * disstep;

    A = cos(theta) * sin(phi);
    B = sin(theta) * sin(phi);
    C = cos(phi);
    D = -d - (A * ptCenter.x + B * ptCenter.y + C * ptCenter.z);
    std::cout << "Plane score " << maxcount << " : ABC : " << A << " , " << B << " , " << C << " , " << D << " theta "
              << theta << " phi " << phi << " d " << d << std::endl;
    foundPlane = Plane(A, B, C, D, ptCenter);
    foundPlane.theta = theta;
    foundPlane.phi = phi;
    foundPlane.dis = d;
    foundPlane.count = maxcount;
    //释放cube
    for (int i = 0; i < thetas; ++i) {
        int **row = cube[i];
        for (int j = 0; j < phis; ++j) {
            int *col = row[j];
            delete[] col;
        }
        delete[] row;
    }
    delete[] cube;
}

int main() {

//    ///create window
//    pangolin::CreateWindowAndBind("Main", 640, 480);
//    glEnable(GL_DEPTH_TEST);
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    pangolin::OpenGlRenderState s_cam(
//            pangolin::ProjectionMatrix(640, 480, 420, 420, 320, 320, 0.2, 100),
//            pangolin::ModelViewLookAt(2, 0, 2, 0, 0, 0, pangolin::AxisY)
//    );
//    pangolin::Handler3D handler(s_cam);
//    pangolin::View &d_cam = pangolin::CreateDisplay()
//            .SetBounds(0.0, 1.0, 0.0, 1.0, -640.0f / 480.0f)
//            .SetHandler(&handler);

    ///Load data
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//    // Generate pointcloud data
//    cloud->width = 100000;
//    cloud->height = 1;
//    for(size_t i=0; i<cloud->points.size();i++)
//    {
//        cloud->points[i].x = 1024.0f * rand()/(RAND_MAX+1.0f);
//        cloud->points[i].y = 1024.0f * rand()/(RAND_MAX+1.0f);
//        cloud->points[i].z = 1024.0f * rand()/(RAND_MAX+1.0f);
//    }
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
            vector<Point> thisCluster;
            for (int index = 0; index < clusters[ci].indices.size(); index++) {
                thisCluster.push_back(Point(cloud->points[clusters[ci].indices[index]].x,
                                            cloud->points[clusters[ci].indices[index]].y,
                                            cloud->points[clusters[ci].indices[index]].z));
            }
            double A, B, C, D;
            Plane thisPlane(Point(0, 0, 0));
            double startTime = clock();
            HoughTransform(thisCluster, A, B, C, D, thisPlane,1);
            Plane thisPlane2(Point(0, 0, 0));
            HoughTransform(thisCluster, A, B, C, D, thisPlane2,2);
            double endTime = clock();
            double timeUsed = double(endTime - startTime) / CLOCKS_PER_SEC;
            cout<<"Hough Transform x 2 time consumed : "<<timeUsed<<endl;
            if(thisPlane.count>thisPlane2.count)
                planes.push_back(thisPlane);
                else
                planes.push_back(thisPlane2);
            ///no limit case
            //startTime = clock();
//            Plane thisPlane3(Point(0, 0, 0));
//            HoughTransform(thisCluster, A, B, C, D, thisPlane3, 3);
//            endTime = clock();
//            timeUsed = double(endTime - startTime) / CLOCKS_PER_SEC;
//            cout << "Hough Transform Full time cosumed : " << timeUsed << endl;
//            cout << "plane 1 contains " << thisPlane.count << " plane 2 " << thisPlane2.count << " plane 3 "
//                 << thisPlane3.count << endl;

            //planes.push_back(thisPlane3);
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

            //需要绘制的东西写在这里
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
                    glEnd();//点设置的结束
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
