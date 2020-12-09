#include <iostream>
//#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstring>
#include <pangolin/pangolin.h>

using namespace std;
#define PI 3.141592653

class Point {
public:
    double x, y, z;
    int planID;
    double distance;

    Point(double ix, double iy, double iz) :
            x(ix), y(iy), z(iz) { planID = -1; distance = -1; }

    Point operator-(const Point &pt) const {
        return Point(x - pt.x, y - pt.y, z - pt.z);
    }
};

class Plane {
public:
    double A, B, C, D;
    double phi, theta, dis;
    vector<int> pointList;
    int mergeID;
    Plane(double Ain, double Bin, double Cin, double Din) : A(Ain), B(Bin), C(Cin), D(Din) {mergeID=-1;};
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

    boundingbox() : x_min(std::numeric_limits<double>::max()),
                    y_min(std::numeric_limits<double>::max()),
                    z_min(std::numeric_limits<double>::max()),
                    x_max(-std::numeric_limits<double>::max()),
                    y_max(-std::numeric_limits<double>::max()),
                    z_max(-std::numeric_limits<double>::max()) {}

    Point center() const {
        return Point((x_max + x_min) / 2.0, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0);
    }
};

void calcBoundBox(vector<Point> points, boundingbox &box) {
    for (int i = 0, n = points.size(); i < n; i++) {
        auto thisP = points[i];
        if (thisP.x < box.x_min)
            box.x_min = thisP.x;
        if (thisP.x > box.x_max)
            box.x_max = thisP.x;
        if (thisP.y < box.y_min)
            box.y_min = thisP.y;
        if (thisP.y > box.y_max)
            box.y_max = thisP.y;
        if (thisP.z < box.z_min)
            box.z_min = thisP.z;
        if (thisP.z > box.z_max)
            box.z_max = thisP.z;
    }
}

void HoughTransfrom(std::vector<Point> &inputs, vector<Plane>& planeList, double &maxA, double &maxB, double &maxC, double &maxD) {
    int n = inputs.size();
    if (n < 3) {
        cout << "point num less than 3" << endl;
        return;
    }

    double theta_start = 0, theta_end = PI;
    double phi_start = 0, phi_end = PI;
    //double phi_start = 0.25*PI, phi_end = 0.75*PI;
    double anglestep = PI / 90, disstep = 0.01;

    boundingbox box;
    calcBoundBox(inputs, box);
    double d_start = -box.diag() / 2.0, d_end = box.diag() / 2.0;

    int thetas = ceil((theta_end - theta_start) / anglestep);
    int phis = ceil((phi_end - phi_start) / anglestep);
    int dises = ceil(box.diag() / disstep);

    int ***cube = new int **[thetas];
    for (int i = 0; i < thetas; i++) {
        cube[i] = new int *[phis];
        for (int j = 0; j < phis; j++) {
            cube[i][j] = new int[dises];
            memset(cube[i][j], 0, sizeof(int) * dises);
        }
    }

    //cos(theta)sin(phi)X+sin(theta)sin(phi)Y+cos(phi)Z = D
    Point ptCenter = box.center();
    for (int i = 0; i < n; i++) {
        if (i % 5000 == 0)
            cout << "processing point " << i << endl;
        const Point &ptOrigin = inputs[i];
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
                ++(col[d_index]);

                phi += anglestep;
                if (phi > phi_end)
                    break;
            }
            theta += anglestep;
            if (phi > theta_end)
                break;
        }
    }

    int buf = 1;
    int maxcount = 0;
    int imax, jmax, kmax;
    vector<vector<double>> selectedPlanes;
    double countThreshold = 10000;
    for (int i = 0; i < thetas; ++i) {
        for (int j = 0; j < phis; ++j) {
            for (int k = buf; k < dises-buf; ++k)
            {
                int count = 0;
                ////wrap nearby cubes i-buf, i, i+buf
                for (int x = i - buf; x <= i + buf; ++x)
                    for (int y = j - buf; y <= j + buf; ++y)
                        for (int z = k - buf; z <= k + buf; ++z) {
                            count += cube[x < 0 ? x + thetas : x % thetas][y < 0 ? y + phis : y % phis][z];
                        }
                //count = cube[i][j][k];
                if (count > maxcount) {
                    imax = i;
                    jmax = j;
                    kmax = k;
                    maxcount = count;
                }
//                if (count > countThreshold) {
//                    vector<double> thisPlane;
//                    thisPlane.push_back(i);
//                    thisPlane.push_back(j);
//                    thisPlane.push_back(k);
//                    thisPlane.push_back(count);
//                    selectedPlanes.push_back(thisPlane);
//                }
            }
        }
    }
    double theta = theta_start + imax * anglestep;
    double phi = phi_start + jmax * anglestep;
    double d = d_start + kmax * disstep;

    double A = cos(theta) * sin(phi);
    double B = sin(theta) * sin(phi);
    double C = cos(phi);
    double D = -d - (A * ptCenter.x + B * ptCenter.y + C * ptCenter.z);
    std::cout <<"max plane"<< A << " , " << B << " , " << C << " , " << D << std::endl;
    maxA = A; maxB = B; maxC = C; maxD = D;

//    for (int pi = 0; pi < selectedPlanes.size(); pi++) {
//        double thisTheta = theta_start + selectedPlanes[pi][0] * anglestep;
//        double thisPhi = phi_start + selectedPlanes[pi][1] * anglestep;
//        double thisd = d_start + selectedPlanes[pi][2] * disstep;
//        double thisA = cos(thisTheta) * sin(thisPhi);
//        double thisB = sin(thisTheta) * sin(thisPhi);
//        double thisC = cos(thisPhi);
//        double thisD = -thisd - (thisA * ptCenter.x + thisB * ptCenter.y + thisC * ptCenter.z);
//        //cout<<"Plane contains "<<selectedPlanes[pi][3]<<" points --- "<<thisA<<" "<<thisB<<" "<<thisC<<" "<<thisD<<endl;
//        Plane newPlane = Plane(thisA,thisB,thisC,thisD);
//        newPlane.phi = thisPhi; newPlane.theta = thisTheta; newPlane.dis = thisd;
//        planeList.push_back(newPlane);
//    }
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
    std::cout << "Hello, World!" << std::endl;
    //Load data
    std::vector<Point> laserPoints;
//    for(int i = 0; i<1000000; i++)
//    {
//        Point thisP = Point(0,0,0);
//        laserPoints.push_back(thisP);
//    }
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
    for (int i = 0; i < num; i++) {
//        if(laserPoints.size()>=20000)
//            break;
        Point newP = Point(*px, *py, *pz);
        laserPoints.push_back(newP);
        px += 4;
        py += 4;
        pz += 4;
        pr += 4;
    }
    fclose(fstream);
    cout << "loaded " << num << " points" << endl;
    cout << "loaded " << laserPoints.size() << " points" << endl;

    //run plan finding
    vector<Plane> foundPlanes;
    double maxA,maxB,maxC,maxD;
    HoughTransfrom(laserPoints, foundPlanes, maxA, maxB, maxC, maxD);
    cout<<"found "<<foundPlanes.size()<<" plane"<<endl;

    ///plane corner points
    vector<Point> corners;
    double cx = -2.5, cy =-1.5;
    double cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));
    cx = -2.5; cy = 1.5; cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));
    cx = 2.5; cy = -1.5; cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));
    cx = 2.5; cy = 1.5; cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));

    //TODO to fix
//    ///merge planes
//    vector<Plane> validPlane;
//    for (int pli = 0; pli < foundPlanes.size(); pli++) {
//        if (foundPlanes[pli].mergeID == -1) {
//            for (int plj = pli + 1; plj < foundPlanes.size(); plj++) {
//                double theta_dif = foundPlanes[pli].theta - foundPlanes[plj].theta;
//                double phi_dif = foundPlanes[pli].phi - foundPlanes[plj].phi;
//                double dis_dif = foundPlanes[pli].dis - foundPlanes[plj].dis;
//                cout << pli << " " << plj << " // ";
//                cout << "theta " << foundPlanes[pli].theta << " - " << foundPlanes[plj].theta << " = "
//                     << theta_dif
//                     << " | " << "phi " << foundPlanes[pli].phi << " - " << foundPlanes[plj].phi << " = "
//                     << phi_dif
//                     << " | " << "dis " << foundPlanes[pli].dis << " - " << foundPlanes[plj].dis << " = "
//                     << dis_dif;
//                if (abs(theta_dif) <= 0.1 && abs(phi_dif) <= 0.1 && abs(dis_dif) <= 1) {
//                    foundPlanes[plj].mergeID = pli;
//                    cout<<" !!! Plane "<<pli<<" merge plane "<<plj;
//                    for (int pi = 0; pi < foundPlanes[plj].pointList.size(); pi++) {
//                        foundPlanes[pli].pointList.push_back(foundPlanes[plj].pointList[pi]);
//                    }
//                }
//                cout<<endl;
//            }
//        }
//    }
//    int validPlaneCount = 0;
//    for (int pli = 0; pli < foundPlanes.size(); pli++)
//    {
//        if(foundPlanes[pli].mergeID == -1)
//        {
//            validPlane.push_back(foundPlanes[pli]);
//            validPlaneCount ++;
//            cout<<"ID "<<pli<<" "<<foundPlanes[pli].theta<<" "<<foundPlanes[pli].phi<<" "<<foundPlanes[pli].dis<<" | "<<foundPlanes[pli].A<<" "<<foundPlanes[pli].B<<" "<<foundPlanes[pli].C<<" "<<foundPlanes[pli].D<<endl;
//        }
//    }
//
//    ///assign points to plane
//    for (int pi = 0; pi < laserPoints.size(); pi++) {
//        double distance = 65535;
//        cout<<endl;
//        for (int pli = 0; pli < validPlane.size(); pli++) {
//            double fenmu = sqrt(validPlane[pli].A * validPlane[pli].A + validPlane[pli].B * validPlane[pli].B + validPlane[pli].C * validPlane[pli].C);
//            double distance_tmp =
//                    abs(validPlane[pli].A * laserPoints[pi].x + validPlane[pli].B * laserPoints[pi].y + validPlane[pli].C * laserPoints[pi].z + validPlane[pli].D) / fenmu;
//            //cout<<distance_tmp<<" ";
//            if (distance_tmp < 0.10) {
//                int Test = 0;
//                if (distance_tmp < distance) {
//                    distance = distance_tmp;
//                    laserPoints[pi].planID = pli;
//                }
//            }
//        }
//        //cout<<"point assigned to "<<laserPoints[pi].planID<<endl;
//    }
    double down = sqrt(maxA * maxA + maxB * maxB + maxC * maxC);
    for (int pi = 0; pi < laserPoints.size(); pi++) {
        double up = abs(maxA * laserPoints[pi].x + maxB * laserPoints[pi].y + maxC * laserPoints[pi].z + maxD);
        double distance = up / down;
        if (distance < 0.05) {
            laserPoints[pi].planID = 1;
        }
        cout << "point " << laserPoints[pi].x << " " << laserPoints[pi].y << " " << laserPoints[pi].z << " dis "
             << distance << "| plain ID "<<laserPoints[pi].planID<<endl;
    }
    cout<<"max Plane "<<maxA<<" "<<maxB<<" "<<maxC<<" "<<maxD<<endl;

    //show
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

        ///Draw
        for (int pi = 0; pi < laserPoints.size(); pi++) {
            glBegin(GL_POINTS);
            switch (laserPoints[pi].planID){
                case -1 :
                    glColor3f(0.8, 0.8, 0.8);
                    break;
                case 1 :
                    glColor3f(1, 0, 0);
                    break;
                case 2 :
                    glColor3f(0, 1, 0);
                    break;
                case 3 :
                    glColor3f(0, 0, 1);
                    break;
                case 4 :
                    glColor3f(1, 1, 0);
                    break;
            }
            glVertex3f(laserPoints[pi].x, laserPoints[pi].y, laserPoints[pi].z);
            glEnd();//点设置的结束
        }
        //lines
        glBegin(GL_LINES);
        glLineWidth(2.0);
        glColor3f(0,0,1);
        glVertex3d(corners[0].x,corners[0].y,corners[0].z);
        glVertex3d(corners[1].x,corners[1].y,corners[1].z);
        glVertex3d(corners[0].x,corners[0].y,corners[0].z);
        glVertex3d(corners[2].x,corners[2].y,corners[2].z);
        glVertex3d(corners[1].x,corners[1].y,corners[1].z);
        glVertex3d(corners[3].x,corners[3].y,corners[3].z);
        glVertex3d(corners[2].x,corners[2].y,corners[2].z);
        glVertex3d(corners[3].x,corners[3].y,corners[3].z);
        glEnd();
        ///End of Draw



        pangolin::FinishFrame();
    }

    return 0;
}
