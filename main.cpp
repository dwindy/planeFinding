#include <iostream>
//#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstring>
#include <pangolin/pangolin.h>

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
typedef Point Vector;

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

void HoughTransform(const std::vector<Point>& input, double& A, double& B, double& C, double& D)
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
                }
            }

    double theta = theta_start + xmax*anglestep;
    double phi = phi_start + ymax*anglestep;
    double d = d_start + zmax*disstep;

    A = cos(theta)*sin(phi);
    B = sin(theta)*sin(phi);
    C = cos(phi);
    D = -d - (A*ptCenter.x + B*ptCenter.y+C*ptCenter.z);
    std::cout <<"Plane"<< A << " , " << B << " , " << C << " , "<< D << std::endl;
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
        Point newP = Point((double)*px, (double)*py, (double)*pz);
        if(*px>0 && *px<25 && abs(*py)<6 )//&& *pz<=-1.4)
            laserPoints.push_back(newP);
        px += 4;
        py += 4;
        pz += 4;
        pr += 4;
    }
    fclose(fstream);
    cout << "loaded " << num << " points" << endl;
    cout << "loaded " << laserPoints.size() << " points" << endl;

///write points
//    std::ofstream  outfile("new.txt",std::ofstream::binary);
//    int writeCount = 0;
//    for(int i = 0; i<laserPoints.size();i++)
//    {
//        if(abs(laserPoints[i].x)<10 && abs(laserPoints[i].y)<5 && abs(laserPoints[i].z<3))
//        {
//            outfile<<laserPoints[i].x<<" "<<laserPoints[i].y<<" "<<laserPoints[i].z<<endl;
//            writeCount++;
//        }
//    }
//    outfile.close();
//    cout<<"wrote "<<writeCount<<" points"<<endl;
/////fake data
//    std::ofstream  outfile("fakedata.txt",std::ofstream::binary);
//    int writeCount = 0;
//    for (float i = -0.25; i < 0.25; i = i + 0.01)
//        for (float j = -0.15; j < 0.15; j = j + 0.01) {
//            double y = (rand()%10+1)/1000.0;
//            cout<<y<<endl;
//            outfile << i << " " << j << " " << y << endl;
//        }
////    for (float j = -0.05; j < 0.05; j = j + 0.01) {
////        for (float k = -0.1; k < 0.1; k = k + 0.01) {
////            outfile << 0 << " " << j << " " << k << endl;
////        }
////    }

/////load from txt
//    std::ifstream infile;
//    infile.open("fakedata.txt");
//    vector<Point> laserPoints;
//    if(!infile)
//    {
//        cout<<"read error"<<endl;
//        return 0;
//    }
//    double x,y,z;
//    while(infile>>x>>y>>z)
//    {
//        Point newP = Point(x,y,z);
//        laserPoints.push_back(newP);
//        //cout<<"read "<<newP.x<<" "<<newP.y<<" "<<newP.z<<endl;
//    }

    ///hough plane
    double maxA,maxB,maxC,maxD;
    HoughTransform(laserPoints,maxA,maxB,maxC,maxD);

    ///plane corner points
    vector<Point> corners;
    double cx = -0.25, cy =-0.15;
    double cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));
    cx = -0.25; cy = 0.15; cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));
    cx = 0.25; cy = -0.15; cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));
    cx = 0.25; cy = 0.15; cz = (maxA*cx + maxB*cy + maxD)/-maxC;
    corners.push_back(Point(cx,cy,cz));

    ///distance from point to plane
    double down = sqrt(maxA*maxA+maxB*maxB+maxC*maxC);
    for(int pi=0;pi<laserPoints.size();pi++)
    {
        double up = abs(maxA * laserPoints[pi].x + maxB * laserPoints[pi].y + maxC*laserPoints[pi].z + maxD);
        double distance = up / down;
        if(distance < 0.10)
            laserPoints[pi].planID=1;
        cout<<"point "<<laserPoints[pi].x<<" "<<laserPoints[pi].y<<" "<<laserPoints[pi].z<<" dis "<<distance<<" | planeID "<<laserPoints[pi].planID<<endl;
    }
    cout<<"max Plane "<<maxA<<" "<<maxB<<" "<<maxC<<" "<<maxD<<endl;

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
            glVertex3d(laserPoints[pi].x, laserPoints[pi].y, laserPoints[pi].z);
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
        //

        pangolin::FinishFrame();
    }

    return 0;
}
