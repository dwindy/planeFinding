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

using namespace std;
#define PI 3.141592653

class Point
{
public:
    double x, y, z;
    int planID;
    Point(double ix,double iy,double iz) :
            x(ix), y(iy), z(iz){planID=-1;}

    Point operator-(const Point& pt) const
    {
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
    Plane(double Ain, double Bin, double Cin, double Din) : A(Ain), B(Bin), C(Cin), D(Din) {mergeID=-1;count=-1;};
};

class boundingbox
{
public:
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;

public:
    double diag() const
    {
        double dx = x_max - x_min;
        double dy = y_max - y_min;
        double dz = z_max - z_min;

        return sqrt(dx*dx + dy*dy + dz*dz);
    }

    boundingbox():
            x_min(std::numeric_limits<double>::max()),
            y_min(std::numeric_limits<double>::max()),
            z_min(std::numeric_limits<double>::max()),
            x_max(-std::numeric_limits<double>::max()),
            y_max(-std::numeric_limits<double>::max()),
            z_max(-std::numeric_limits<double>::max())
    {}

    Point center() const
    {
        return Point((x_max + x_min) / 2.0,(y_min+y_max) / 2.0, (z_min+z_max) / 2.0);
    }
};

void calcboundbox(const std::vector<Point>& input, boundingbox& box)
{
    for (int i = 0, n = input.size(); i < n;++i)
    {
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

void HoughTransform(const std::vector<Point>& input, vector<Plane>& foundPlanes, double& A, double& B, double& C, double& D)
{
    int n = input.size();
    if (n < 3)
        return;

    double theta_start=0, theta_end=PI;
    double phi_start=0, phi_end=PI;
    //double phi_start = 0.25*PI, phi_end = 0.75*PI;
    double anglestep=PI/90, disstep=0.05;

    boundingbox box;
    calcboundbox(input, box);
    double d_start = -box.diag() / 2.0, d_end = box.diag() / 2.0;

    int thetas = ceil((theta_end - theta_start) / anglestep);
    int phis = ceil((phi_end - phi_start) / anglestep);
    int dises = ceil( box.diag()/disstep);

    int*** cube = new int**[thetas];
    for (int i = 0; i < thetas;++i)
    {
        cube[i] = new int*[phis];
        for (int j = 0; j < phis; ++j)
        {
            cube[i][j] = new int[dises];
            memset(cube[i][j], 0, sizeof(int)*dises);
        }
    }

    //cos(theta)sin(phi)X+sin(theta)sin(phi)Y+cos(phi)Z = D
    Point ptCenter = box.center();
    for (int i = 0; i < n;++i)
    {
        if(i%5000==0)
            cout<<"dealing with point "<<i<<endl;
        const Point& ptOrigin = input[i];
        Point point = ptOrigin - ptCenter;

        double theta = theta_start;
        for(int j = 0; j < thetas; ++j)
        {
            int** row = cube[j];
            double phi = phi_start;
            for (int k = 0; k < phis; ++k)
            {
                int* col = row[k];

                double sinphi = sin(phi);
                double d = cos(theta)*sinphi*point.x + sin(theta)*sinphi*point.y + cos(phi)*point.z;

                int d_index = floor((d - d_start) / disstep);
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
    for (int i = 0; i < thetas;++i)
        for (int j = 0; j < phis; ++j)
            for (int k = buf; k < dises - buf;++k)
            {
                int count = 0;
                for (int x = i - buf; x <= i + buf; ++x)
                    for (int y = j - buf; y <= j + buf; ++y)
                        for (int z = k - buf; z <= k + buf; ++z)
                        {
                            count += cube[x<0?x+thetas:x%thetas][y<0?y+phis:y%phis][z];
                        }
                if (count > maxcount)
                {
                    xmax = i;
                    ymax = j;
                    zmax = k;
                    maxcount = count;
                    //cout<<"new max count "<<maxcount<<endl;
                }
                if(count > countThres)
                {
                    double thisTheta = theta_start + i * anglestep;
                    double thisPhi = phi_start + j * anglestep;
                    double thisDis = d_start + k * disstep;
                    int thisCount = count;
                    double thisA = cos(thisTheta)*sin(thisPhi);
                    double thisB = sin(thisTheta)*sin(thisPhi);
                    double thisC = cos(thisPhi);
                    double thisD = -thisDis - (thisA*ptCenter.x + thisB*ptCenter.y + thisC*ptCenter.z);
                    Plane thisPlane = Plane(thisA, thisB, thisC, thisD);
                    thisPlane.theta = thisTheta;
                    thisPlane.phi = thisPhi;
                    thisPlane.dis = thisDis;
                    thisPlane.count = thisCount;
                    foundPlanes.push_back(thisPlane);
                }
            }

    double theta = theta_start + xmax*anglestep;
    double phi = phi_start + ymax*anglestep;
    double d = d_start + zmax*disstep;

    A = cos(theta)*sin(phi);
    B = sin(theta)*sin(phi);
    C = cos(phi);
    D = -d - (A*ptCenter.x + B*ptCenter.y+C*ptCenter.z);
    std::cout <<"Plane score "<<maxcount<<" : ABC : "<< A << " , " << B << " , " << C << " , "<< D << " theta "<<theta<<" phi "<<phi<<" d "<<d<<std::endl;
    //释放cube
    for (int i = 0; i < thetas; ++i)
    {
        int** row = cube[i];
        for (int j = 0; j < phis;++j)
        {
            int* col = row[j];
            delete[] col;
        }
        delete[] row;
    }
    delete[] cube;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    //Load data from bin
    std::vector<Point> laserPoints;
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
        Point newP = Point((double) *px, (double) *py, (double) *pz);
        if (*px > 0 && *px < 25 && abs(*py) < 6)
            laserPoints.push_back(newP);
        px += 4;
        py += 4;
        pz += 4;
        pr += 4;
    }
    fclose(fstream);
    cout << "loaded " << num << " points" << endl;
    cout << "loaded " << laserPoints.size() << " points" << endl;

/////generate fake points
//vector<Point> laserPoints;
//    for (float y = -0.25; y < 0.25; y = y + 0.01) {
//        for (float z = -0.15; z < 0.15; z = z + 0.01) {
//            double x = 0.5 + (rand() % 10 + 1) / 1000.0;
//            laserPoints.push_back(Point(x, y, z));
//        }
//    }
//    for (float x = -0.15+0.5; x < 0.15+0.5; x = x + 0.01) {
//        for (float y = -0.25; y < 0.25; y = y + 0.01) {
//            double z = (rand() % 10 + 1) / 1000.0;
//            laserPoints.push_back(Point(x, y, z));
//        }
//    }

    clock_t startTime = clock();
    ///hough plane
    double maxA, maxB, maxC, maxD;
    vector<Plane> foundPlanes;
    HoughTransform(laserPoints, foundPlanes, maxA, maxB, maxC, maxD);
    clock_t endTime = clock();
    double timeUsed = double(endTime - startTime) / CLOCKS_PER_SEC;
    cout << "houghTrans Plane time used " << timeUsed << " sec" << endl;

    ///merge planes
    cout << "found " << foundPlanes.size() << " planes" << endl;
    for (int pli = 0; pli < foundPlanes.size(); pli++)
    {
        Plane thisPlane = foundPlanes[pli];
        if(thisPlane.mergeID==-1)
        {
            //cout<<"plane "<<pli<<" | "<<thisPlane.theta<<" "<<thisPlane.phi<<" "<<thisPlane.dis<<" -> "<<endl;
            for(int plj = pli+1; plj<foundPlanes.size();plj++)
            {
                Plane targetPlane = foundPlanes[plj];
                double thetaDif = abs(thisPlane.theta - targetPlane.theta);
                double phiDif = abs(thisPlane.phi - targetPlane.phi);
                double disDif = abs(thisPlane.dis - targetPlane.dis);
                if (thetaDif <= 0.3 && phiDif <= 0.3 && disDif <= 0.3) {
                    //cout << "target plane "<<plj<<" | "<<targetPlane.theta<<" "<<targetPlane.phi<<" "<<targetPlane.dis<<" merge "<<endl;
                    foundPlanes[plj].mergeID = pli;
                }
            }
        }
    }
    cout<<"selected "<<" plane :"<<endl;
    vector<Plane> selectedPlanes;
    for (int pli = 0; pli < foundPlanes.size(); pli++) {
        if (foundPlanes[pli].mergeID == -1) {
            cout << "plane " << pli << " | " << foundPlanes[pli].A<<" "<< foundPlanes[pli].B<<" "<< foundPlanes[pli].C<<" "<< foundPlanes[pli].D<<" | "<< foundPlanes[pli].theta << " " << foundPlanes[pli].phi << " "
                 << foundPlanes[pli].dis << " -> " << endl;
            selectedPlanes.push_back(foundPlanes[pli]);
        }
    }

///basic plane points
    vector<Point> uniplanecorners;
    double unitA = 0, unitB = 0, unitC = 1, unitD = 0;
    double cx = -0.03, cy =-0.05;
    double cz = (unitA*cx + unitB*cy + unitD)/-unitC;
    uniplanecorners.push_back(Point(cx,cy,cz));
    cx = -0.03; cy = 0.05; cz = (unitA*cx + unitB*cy + unitD)/-unitC;
    uniplanecorners.push_back(Point(cx,cy,cz));
    cx = 0.03; cy = -0.05; cz = (unitA*cx + unitB*cy + unitD)/-unitC;
    uniplanecorners.push_back(Point(cx,cy,cz));
    cx = 0.03; cy = 0.05; cz = (unitA*cx + unitB*cy + unitD)/-unitC;
    uniplanecorners.push_back(Point(cx,cy,cz));

    ///assign points to planes
    //double down = sqrt(maxA*maxA+maxB*maxB+maxC*maxC);
    for(int pi=0;pi<laserPoints.size();pi++)
    {
        double minDistance = 999;
        for(int pli = 0;pli<selectedPlanes.size();pli++)
        {
            Plane thisPlane = selectedPlanes[pli];
            double up = abs(thisPlane.A * laserPoints[pi].x + thisPlane.B * laserPoints[pi].y + thisPlane.C*laserPoints[pi].z + thisPlane.D);
            double down = sqrt(thisPlane.A*thisPlane.A+thisPlane.B*thisPlane.B+thisPlane.C*thisPlane.C);
            double distance = up/down;
            if(distance < 0.10 && distance<minDistance)
            {
                minDistance = distance;
                laserPoints[pi].planID = pli;
            }
        }
//        double up = abs(maxA * laserPoints[pi].x + maxB * laserPoints[pi].y + maxC*laserPoints[pi].z + maxD);
//        double distance = up / down;
//        if(distance < 0.10)
//            laserPoints[pi].planID=1;
        //cout<<"point "<<laserPoints[pi].x<<" "<<laserPoints[pi].y<<" "<<laserPoints[pi].z<<" dis "<<distance<<" | planeID "<<laserPoints[pi].planID<<endl;
    }
    //cout<<"max Plane "<<maxA<<" "<<maxB<<" "<<maxC<<" "<<maxD<<endl;

    ///maxplane normal vector
    bool flag = false;
    Point selPoint = Point(0,0,0);
    while(!flag)
    {
        int index = rand()%laserPoints.size();
        if(laserPoints[index].planID==1)
        {
            flag = true;
            selPoint.x = laserPoints[index].x;
            selPoint.y = laserPoints[index].y;
            selPoint.z = laserPoints[index].z;
        }
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
        for (int pi = 0; pi < laserPoints.size(); pi++) {
            glBegin(GL_POINTS);
            switch (laserPoints[pi].planID){
                case -1 :
                    glColor3f(0.8, 0.8, 0.8);
                    break;
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
                    glColor3f(0, 1, 1);
                    break;
            }
            glVertex3d(laserPoints[pi].x, laserPoints[pi].y, laserPoints[pi].z);
            glEnd();//点设置的结束
        }
        //box of max plane
        glBegin(GL_LINES);
        glLineWidth(2.0);
        glColor3f(1,1,1);
        glVertex3d(selPoint.x,selPoint.y,selPoint.z);
        glVertex3d(selPoint.x+maxA,selPoint.y+maxB,selPoint.z+maxC);
        glEnd();
        //coordinate frame
        glBegin(GL_LINES);
        glLineWidth(5.0);
        glColor3f(1,0,1);
        glVertex3d(0.1,0,0);
        glVertex3d(0,0,0);
        glColor3f(0,1,0);
        glVertex3d(0,0.1,0);
        glVertex3d(0,0,0);
        glColor3f(0,0,1);
        glVertex3d(0,0,0.1);
        glVertex3d(0,0,0);
        glEnd();
        //box of unit plane
        glBegin(GL_LINES);
        glLineWidth(10.0);
        glColor3f(1,1,1);
        glVertex3d(uniplanecorners[0].x,uniplanecorners[0].y,uniplanecorners[0].z);
        glVertex3d(uniplanecorners[1].x,uniplanecorners[1].y,uniplanecorners[1].z);
        glVertex3d(uniplanecorners[0].x,uniplanecorners[0].y,uniplanecorners[0].z);
        glVertex3d(uniplanecorners[2].x,uniplanecorners[2].y,uniplanecorners[2].z);
        glVertex3d(uniplanecorners[1].x,uniplanecorners[1].y,uniplanecorners[1].z);
        glVertex3d(uniplanecorners[3].x,uniplanecorners[3].y,uniplanecorners[3].z);
        glVertex3d(uniplanecorners[2].x,uniplanecorners[2].y,uniplanecorners[2].z);
        glVertex3d(uniplanecorners[3].x,uniplanecorners[3].y,uniplanecorners[3].z);
        glColor3f(1,1,1);
        glVertex3d(0,0,0);
        glVertex3d(0,0,0.05);
        glEnd();
        //

        pangolin::FinishFrame();
    }

    return 0;
}
