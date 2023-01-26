#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <lbfgs.h>
#include <GL/glut.h>
#include <vector>
#include <random>
using namespace std;
#define SEGMENT 10 // 一本の糸のセグメント数
#define YARN 20 //糸の総数
#define N   SEGMENT*YARN*2 //全要素（高さと傾き）
#define DIMENSION 10 //分割数
#define SPLINE SEGMENT*YARN*DIMENSION //スプライン曲線全要素
#define RADIUS 0.5 //糸の半径
#define INTERVAL RADIUS*2 //(仮) //同じ方向の糸同士の間隔
#define Ks 2 //
#define Kc 200 //
#define CIRCLE 30 //回転の分割数
#define RATIO INTERVAL*YARN/(2*SEGMENT)  //比率
#define PLY 3 //より糸の本数
#define Rply 0.27 //より糸と糸の距離
#define Aply 5.0 //より糸の螺旋のピッチ
#define Tply 0 //より糸の螺旋の初期位相
#define Rradius 1.0//より糸の半径
#define FIBER 200 //繊維の本数
#define Fradius 0.05 //繊維の半径
#define Afiber 5.0 //繊維の螺旋のピッチ
#define Rayply 360 //より糸のレイ
#define length 10 //より糸の長さ
#define plypitch 20 //より糸の螺旋ピッチ
vector<GLfloat> K(2);
static const double pi = 3.141592653589793;
lbfgsfloatval_t radius = (lbfgsfloatval_t)RADIUS;
lbfgsfloatval_t circle = (lbfgsfloatval_t)CIRCLE;
vector<vector<GLfloat>> F(FIBER , vector<GLfloat>(2));

vector<vector<GLfloat>> fibersample(FIBER,vector<GLfloat>(3));
vector<vector<vector<GLfloat>>> fibersample2(CIRCLE,vector<vector<GLfloat>>(FIBER,vector<GLfloat>(3)));
vector<vector<GLfloat>> RayPly(Rayply, vector<GLfloat>(3));
vector<vector<GLfloat>> Visible(Rayply, vector<GLfloat>(3));
vector<vector<GLfloat>> normal(Rayply, vector<GLfloat>(3));
vector<vector<GLfloat>> PLYnormal(Rayply, vector<GLfloat>(3));
int Shadow[Rayply];
float RayPly2[Rayply][3]; 
float RayPly10[100][Rayply][3];
float Tangent[Rayply][3][length];


vector<GLfloat> white = { 1.0,1.0,1.0,1.0 };
vector<GLfloat> black = { 0.0,0.0,0.0,1.0 };
vector<GLfloat> half = { 0.5,0.5,0.5,1.0 };
vector<GLfloat> half2 = { 0.1,0.1,0.1,1.0 };
vector<GLfloat> red = { 1.0,0.0,0.0,0.5 };
vector<vector<GLfloat>> red2(PLY,vector<GLfloat>(4));
vector<vector<GLfloat>> red3(FIBER,vector<GLfloat>(4));
vector<GLfloat> blue = { 0.0,0.0,1.0,1.0 };
vector<vector<GLfloat>> blue2(PLY, vector<GLfloat>(4));
vector<vector<GLfloat>> blue3(FIBER, vector<GLfloat>(4));
vector<GLfloat> green = { 0.0,1.0,0.0,1.0 };
vector<GLfloat> pink = { 0.9,0.43,0.43,1.0 };

vector<GLfloat> light1pos = { 0.0,0.0,15.0,0,1.0 };
vector<GLfloat> lightDiffuse = {0.9,0.9,0.9 };
vector<GLfloat> lightAmbient = { 0.2,0.2,0.2 };
vector<GLfloat> lightSpecular = { 0.01,0.01,0.01 };
vector<GLfloat> spotDirrection = { -0.5,0.0,-1.0 };

vector<vector<vector<GLfloat>>> V1(CIRCLE,vector<vector<GLfloat>>(SPLINE,vector<GLfloat>(3)));
vector<vector<vector<GLfloat>>> V2(CIRCLE, vector<vector<GLfloat>>(SPLINE, vector<GLfloat>(3)));
vector<vector<vector<vector<GLfloat>>>> V3(PLY, vector<vector<vector<GLfloat>>>(CIRCLE, vector<vector<GLfloat>>(SPLINE,vector<GLfloat>(3))));
vector<vector<vector<vector<GLfloat>>>> V4(PLY, vector<vector<vector<GLfloat>>>(CIRCLE, vector<vector<GLfloat>>(SPLINE, vector<GLfloat>(3))));
vector< vector<vector<vector<vector<GLfloat>>>>> V5(PLY, vector<vector<vector<vector<GLfloat>>>>(FIBER, vector<vector<vector<GLfloat>>>(CIRCLE,vector<vector<GLfloat>>(SPLINE,vector<GLfloat>(3)))));
vector< vector<vector<vector<vector<GLfloat>>>>> V6(PLY, vector<vector<vector<vector<GLfloat>>>>(FIBER, vector<vector<vector<GLfloat>>>(CIRCLE, vector<vector<GLfloat>>(SPLINE, vector<GLfloat>(3)))));
float tex1D[Rayply][3];
float tex1D10[10][Rayply][3];
GLfloat r = 0.0;
/* ドラッグ開始位置 */
static int cx, cy;

/* マウスの絶対位置→ウィンドウ内での相対位置の換算係数 */
static double sx, sy;

/* マウスの相対位置→回転角の換算係数 */
#define SCALE (2.0 * 3.14159265358979323846)

/* 回転の初期値 (クォータニオン) */
static double cq[4] = { 1.0, 0.0, 0.0, 0.0 };

/* ドラッグ中の回転 (クォータニオン) */
static double tq[4];

/* 回転の変換行列 */
static double rt[16];

GLfloat theta = 0.0;

class vec3 {
private:
public:
    float x, y, z;
};


GLfloat thetaX(GLfloat theta, GLfloat x, GLfloat y, GLfloat z) {
    GLfloat Ri = sqrt(x * x + y * y + z * z);
    GLfloat ti = atan2(y, x);
    GLfloat alhfa = 1.0;

    return Ri * cos(theta + ti);
    
}

GLfloat thetaY(GLfloat theta, GLfloat x, GLfloat y, GLfloat z) {
    GLfloat Ri = sqrt(x * x + y * y + z * z);
    GLfloat ti = atan2(y, x);
    GLfloat alhfa = 1.0;

    return Ri * sin(theta + ti);
     

}

GLfloat thetaZ(GLfloat theta, GLfloat x, GLfloat y, GLfloat z) {
    GLfloat Ri = sqrt(x * x + y * y + z * z);
    GLfloat ti = atan2(y, x);
    GLfloat alhfa = 0.0;
    return alhfa * theta / (2 * pi);
}

GLfloat tangentX(GLfloat theta, GLfloat x, GLfloat y,GLfloat z) {
    GLfloat Ri = sqrt(x * x + y * y + z * z);
    GLfloat ti = atan2(y, x);
    GLfloat alhfa = 1.0;

    return -Ri * sin(theta + ti);

}

GLfloat tangentY(GLfloat theta, GLfloat x, GLfloat y,GLfloat z) {
    GLfloat Ri = sqrt(x * x + y * y + z * z);
    GLfloat ti = atan2(y, x);
    GLfloat alhfa = 1.0;

    return Ri * cos(theta + ti);


}

GLfloat tangentZ(GLfloat theta, GLfloat x, GLfloat y,GLfloat z) {
    GLfloat Ri = sqrt(x * x + y * y + z * z);
    GLfloat ti = atan2(y, x);
    GLfloat alhfa = 0.0;
    return alhfa / (2 * pi);
}

/*
** クォータニオンの積 r <- p x q
*/
void qmul(double r[], const double p[], const double q[])
{
    r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
    r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
    r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
    r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

/*
** 回転の変換行列 r <- クォータニオン q
*/
void qrot(double r[], double q[])
{
    double x2 = q[1] * q[1] * 2.0;
    double y2 = q[2] * q[2] * 2.0;
    double z2 = q[3] * q[3] * 2.0;
    double xy = q[1] * q[2] * 2.0;
    double yz = q[2] * q[3] * 2.0;
    double zx = q[3] * q[1] * 2.0;
    double xw = q[1] * q[0] * 2.0;
    double yw = q[2] * q[0] * 2.0;
    double zw = q[3] * q[0] * 2.0;

    r[0] = 1.0 - y2 - z2;
    r[1] = xy + zw;
    r[2] = zx - yw;
    r[4] = xy - zw;
    r[5] = 1.0 - z2 - x2;
    r[6] = yz + xw;
    r[8] = zx + yw;
    r[9] = yz - xw;
    r[10] = 1.0 - x2 - y2;
    r[3] = r[7] = r[11] = r[12] = r[13] = r[14] = 0.0;
    r[15] = 1.0;
}


void idle(void) {
    glutPostRedisplay();
}

void display(void) {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor4fv(&white[0]);
    glLoadIdentity();



    gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    glLightfv(GL_LIGHT1, GL_POSITION, light1pos.data());



    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, pink.data());
   /*回転*/
    glMultMatrixd(rt);
    //glRotated(r, 0.0, 1.0, 0.0);
    /*
    glBegin(GL_QUADS);
    int oo = 0;
    for (int j = 0; j < 9; j++) {
        glNormal3fv(&tex1D10[j][0][0]);
        glVertex3fv(&RayPly10[j][1 + oo][0]);
        glVertex3fv(&RayPly10[j][0 + oo][0]);
        glVertex3fv(&RayPly10[j + 1][1 + oo][0]);
        glVertex3fv(&RayPly10[j + 1][2 + oo][0]);
        oo++;
    }
    */
    
    glBegin(GL_QUADS);
    int o = 0;
    for (int j = 0; j < length-1; j++) {
        for (int i = 0; i < Rayply; i++) {
            int a = i + o;
            int b = i + 1 + o;
            int c = i + (Rayply / plypitch) + 1 + o; //法線と形状は対応していなければならない
            int d = i + (Rayply / plypitch) + o;
            if (a >= Rayply) {
                a = a - Rayply;
            }
            if (b >= Rayply) {
                b = b - Rayply;
            }
            if (c >= Rayply) {
                c = c - Rayply;
            }
            if (d >= Rayply) {
                d = d - Rayply;
            }
            glNormal3fv(&tex1D10[j][i][0]);
            glVertex3fv(&RayPly10[j][b][0]);
            glVertex3fv(&RayPly10[j][a][0]);
            glVertex3fv(&RayPly10[j + 1][d][0]);
            glVertex3fv(&RayPly10[j + 1][c][0]);
        }
        o = o + (Rayply / plypitch);
    }
    glEnd();

    //*/
    /*
    glBegin(GL_QUADS);
    for (int i = 0; i < Rayply - 1; i++) {
        glNormal3fv(&tex1D[i][0]);
        glVertex3fv(&RayPly2[i + 1][0]);
        glVertex3fv(&RayPly2[i][0]);
        glVertex3fv(&RayPly[i][0]);
        glVertex3fv(&RayPly[i + 1][0]);
    }
    glNormal3fv(&tex1D[359][0]);
    glVertex3fv(&RayPly2[0][0]);
    glVertex3fv(&RayPly2[359][0]);
    glVertex3fv(&RayPly[359][0]);
    glVertex3fv(&RayPly[0][0]);
    */
    //* FIBER断面
    /*

    glBegin(GL_LINES);
    glColor4fv(&blue[0]);
    glVertex3f(-1.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    glEnd();
    glBegin(GL_LINES);
    glColor4fv(&blue[0]);
    glVertex3f(0.0, -1.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glEnd();
    // */

    glBegin(GL_LINES);
    glColor4fv(&blue[0]);
    for (int k = 0; k < Rayply - 1; k++) {
        glVertex3f(RayPly[k][0], RayPly[k][1], RayPly[k][2]);
        glVertex3f(RayPly[k + 1][0], RayPly[k + 1][1], RayPly[k + 1][2]);
    }
    glVertex3f(RayPly[Rayply - 1][0], RayPly[Rayply - 1][1], RayPly[Rayply - 1][2]);
    glVertex3f(RayPly[0][0], RayPly[0][1], RayPly[0][2]);
    glEnd();
    //

    glBegin(GL_POINTS);
    glColor4fv(&red[0]);
    for (int i = 0; i < fibersample.size(); i++) {
        glColor4fv(&red[0]);
        glVertex3f(thetaX(theta, fibersample[i][0], fibersample[i][1], fibersample[i][2]), thetaY(theta, fibersample[i][0], fibersample[i][1], fibersample[i][2]), thetaZ(theta, fibersample[i][0], fibersample[i][1], fibersample[i][2]));
        fibersample[i][0] = thetaX(theta, fibersample[i][0], fibersample[i][1], fibersample[i][2]);
        fibersample[i][1] = thetaY(theta, fibersample[i][0], fibersample[i][1], fibersample[i][2]);
        fibersample[i][2] = thetaZ(theta, fibersample[i][0], fibersample[i][1], fibersample[i][2]);
    }
    glEnd();


    glBegin(GL_LINES);
    glColor4fv(&red[0]);
    for (int i = 0; i < fibersample.size(); i++) {
        for (int h = 0; h < CIRCLE; h++) {
            lbfgsfloatval_t Gsin = sin((2.0 / CIRCLE) * h * pi);
            lbfgsfloatval_t Gcos = cos((2.0 / CIRCLE) * h * pi);
            fibersample2[h][i][0] = fibersample[i][0] + Gsin * Fradius;
            fibersample2[h][i][1] = fibersample[i][1] + Gcos * Fradius;
            fibersample2[h][i][2] = fibersample[i][2];
            //printf("aaaa %f %f %f\n", fibersample2[h][i][0], fibersample2[h][i][1], fibersample2[h][i][2]);
        }
    }
    for (int i = 0; i < fibersample.size(); i++) {
        for (int h = 0; h < CIRCLE - 1; h++) {
            glVertex3fv(&fibersample2[h][i][0]);
            glVertex3fv(&fibersample2[h + 1][i][0]);
        }
        glVertex3fv(&fibersample2[29][i][0]);
        glVertex3fv(&fibersample2[0][i][0]);
    }
    glEnd();



    glBegin(GL_LINES);
    glColor4fv(&green[0]);
    glPointSize(10.0);
    for (int j = 0; j < Rayply - 1; j++) {
        glColor4fv(&green[0]);
        if (Shadow[j] < 60 || Shadow[j + 1] < 60) {
            glColor4fv(&pink[0]);
        }
        glVertex3f(Visible[j][0], Visible[j][1], Visible[j][2]);
        glVertex3f(Visible[j + 1][0], Visible[j + 1][1], Visible[j + 1][2]);
    }
    //glColor4fv(&pink[0]);
    if (Shadow[359] < 60 || Shadow[0] < 60) {
        glColor4fv(&pink[0]);
    }
    glVertex3f(Visible[359][0], Visible[359][1], Visible[359][2]);
    glVertex3f(Visible[0][0], Visible[0][1], Visible[0][2]);
    glEnd();

    //*/



    glutIdleFunc(idle);


    glutSwapBuffers();

    r += 0.1;
    if (r >= 360) r = 0;

    /*
    r += 0.0;
    if (r >= 360) {
        r = 0;
        theta += pi / 10.0;
    }
    //*/
}



void resize(int w, int h)
{
    /* マウスポインタ位置のウィンドウ内の相対的位置への換算用 */
    sx = 1.0 / (double)w;
    sy = 1.0 / (double)h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(100.0, (double)w / (double)h, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void mouse(int button, int state, int x, int y)
{
    switch (button) {
    case GLUT_LEFT_BUTTON:
        switch (state) {
        case GLUT_DOWN:
            /* ドラッグ開始点を記録 */
            cx = x;
            cy = y;
            /* アニメーション開始 */
            glutIdleFunc(idle);
            break;
        case GLUT_UP:
            /* アニメーション終了 */
            glutIdleFunc(0);
            /* 回転の保存 */
            cq[0] = tq[0];
            cq[1] = tq[1];
            cq[2] = tq[2];
            cq[3] = tq[3];
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }
}

void motion(int x, int y)
{
    double dx, dy, a;

    /* マウスポインタの位置のドラッグ開始位置からの変位 */
    dx = (x - cx) * sx;
    dy = (y - cy) * sy;

    /* マウスポインタの位置のドラッグ開始位置からの距離 */
    a = sqrt(dx * dx + dy * dy);

    if (a != 0.0) {
        /* マウスのドラッグに伴う回転のクォータニオン dq を求める */
        double ar = a * SCALE * 0.5;
        double as = sin(ar) / a;
        double dq[4] = { cos(ar), dy * as, dx * as, 0.0 };

        /* 回転の初期値 cq に dq を掛けて回転を合成 */
        qmul(tq, dq, cq);

        /* クォータニオンから回転の変換行列を求める */
        qrot(rt, tq);
    }
}

void init(void)
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    //glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse.data());
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular.data());
    glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient.data());
    //glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.01);
    //glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spotDirrection);
    //glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 30.0);
    //glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 0.0);

    //glLightfv(GL_LIGHT1, GL_DIFFUSE, white);
    //glLightfv(GL_LIGHT1, GL_SPECULAR, half);
    //glLightfv(GL_LIGHT1, GL_AMBIENT, black);
    /* 回転行列の初期化 */
    qrot(rt, cq);
}




lbfgsfloatval_t p(lbfgsfloatval_t x, lbfgsfloatval_t E, lbfgsfloatval_t Beta) {
    return ((1.0 - 2.0 * E) * pow((exp(1) - exp(x)) / (exp(1) - 1.0), Beta) + E);
}

lbfgsfloatval_t sampling(lbfgsfloatval_t E, lbfgsfloatval_t Beta, lbfgsfloatval_t* en,int i) {
    lbfgsfloatval_t z1 = 10.0;
    lbfgsfloatval_t z2 = 10.0;
    int a = 0;
    lbfgsfloatval_t u = 0.0;
    std::random_device rnd;     // 非決定的な乱数生成器を生成
    std::mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    std::uniform_real_distribution<> rand100(-1, 1);        // [0, 99] 範囲の一様乱数
    std::uniform_real_distribution<> rand101(0, 1);
    while (1) {
        if (a >= 1000000) {
            fibersample.resize(i);
            fibersample2.at(0).resize(i);
            return 1;
        }
        // */
        while (1) {
            z1 = rand100(mt);
            z2 = rand100(mt);
            if (sqrt(z1 * z1 + z2 * z2) < 1.0-Fradius) {
                break;
            }
        }
        u = rand101(mt);
        if (p(sqrt(z1 * z1 + z2 * z2), E, Beta) > u) {
            en[0] = z1;
            en[1] = z2;
            
            if (i > 0) {
                for (int j = 0; j < i; j++) {
                    if (sqrt((F[j][0] - en[0]) * (F[j][0] - en[0]) + (F[j][1] - en[1]) * (F[j][1] - en[1])) > 2*Fradius) { //繊維の直径
                        if (j == i-1) {
                            return 0;
                        }
                    }
                    else {
                        a++;
                        break;
                    }
                }
            }
            else {
                return 0;
            }
        }
    }
}

GLfloat Normal(vector<vector<GLfloat>> a, vector<vector<GLfloat>> b,int i,int j,int k) {
    vector<vector<GLfloat>> c(Rayply, vector<GLfloat>(3));
    GLfloat d = sqrt((b[j][0] - a[i][0]) * (b[j][0] - a[i][0]) + (b[j][1] - a[i][1]) * (b[j][1] - a[i][1]) + (b[j][2] - a[i][2]) * (b[j][2] - a[i][2]));
    c[j][k] = (b[j][k] - a[i][k]) / d;
    return c[j][k];
}


GLfloat kousa(GLfloat xr,GLfloat yr,GLfloat xv,GLfloat yv,GLfloat cx,GLfloat cy,GLfloat r,int t) {
    GLfloat alfa = (yr - yv) / ((xr - xv) + 1e-12);
    GLfloat beta = yr - (alfa*xr);
    GLfloat Ak = 1 + alfa*alfa;
    GLfloat Bk = -2 * cx + 2 * alfa * (beta - cy);
    GLfloat Ck = cx * cx + (beta - cy) * (beta - cy) - r * r;
    GLfloat D = Bk * Bk - 4 * Ak * Ck;
    GLfloat s1 = (-Bk + sqrt(D)) / (2 * Ak);
    GLfloat s2 = (-Bk - sqrt(D)) / (2 * Ak);
     
    if (D > 0) 
    
    {
        if ((xr > xv && xv <= s2) || (xr < xv && xv >= s1)) {
            t++;
        }
        else {
            
        }
        //printf("       [%d] xr %f yr %f xv %f yv %f s1 %f s2 %f\n", t, xr, yr, xv, yv, s1, s2);
        //printf("       a %f b %f\n", alfa, beta);
        //printf("       A %f B %f C %f D %f\n", Ak, Bk, Ck, D);
        
    }
    
    return t;
}




/*
GLfloat kousa(GLfloat xr, GLfloat yr, GLfloat xv, GLfloat yv, GLfloat cx, GLfloat cy, GLfloat r, int t) {
    GLfloat alfa = xv - xr;
    GLfloat beta = yv - yr;
    GLfloat x1 = xr - cx;
    GLfloat y1 = yr - cy;
    GLfloat x2 = xv - cx;
    GLfloat y2 = yv - cy;
    GLfloat pt = -(alfa * x1 + beta * y1) / (alfa * alfa + beta * beta);
    GLfloat dt = (alfa * y1 - beta * x1) * (alfa * y1 - beta * x1) / (alfa * alfa + beta * beta);
    if (pt < 0) {
        dt = x1 * x1 + y1 * y1;
    }
    else if(pt > 1 ) {
        dt = x2 * x2 + y2 * y2;
    }

    if (dt > r * r) {

    }
    else {
        t++;
    }

    return t;
}
*/
void Shadowing() {
    for (int i = 0; i < Rayply; i++) {
        GLfloat th = atan2(Visible[i][1] , Visible[i][0]);
        int t = 0;
        int s = 0;
        for (int k = 0; k < Rayply / 2; k++) {
            GLfloat COS = cos(th - (pi / 2) + (k * 2.0 * pi / Rayply));
            GLfloat SIN = sin(th - (pi / 2) + (k * 2.0 * pi / Rayply));
            GLfloat xr = Visible[i][0] + COS * 2.0 * Fradius;
            GLfloat yr = Visible[i][1] + SIN * 2.0 * Fradius;
            for (int j = 0; j < fibersample.size(); j++) {
                t = kousa(xr,yr,Visible[i][0], Visible[i][1], fibersample[j][0], fibersample[j][1], Fradius, t);
            }
            if (t < 1) {
                s++;
            }
            t = 0;
        }
        
        Shadow[i] = s;
    }
}

/*
GLfloat tangent(int i,int j) {
    Tangent[j][0][i] = Rradius * cos(2 * pi * i / length) + RayPly[j][0];
    Tangent[j][1][i] = Rradius * sin(2 * pi * i / length) + RayPly[j][1];
    Tangent[j][2][i] = Afiber * 2 * pi * i / length + RayPly[j][2];
}

GLfloat Tangentmake() {
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < Rayply; j++) {
            tangent(i, j);
        }
    }
}
*/

int main(int argc, char **argv)
{
    //objective_function obj;
    //obj.run(N);
    lbfgsfloatval_t S = 0.0;
    lbfgsfloatval_t s = 1.0 / (lbfgsfloatval_t)DIMENSION;
    lbfgsfloatval_t j = 0.0;
    lbfgsfloatval_t interval = (lbfgsfloatval_t)INTERVAL;
    lbfgsfloatval_t V = 0.0;
    lbfgsfloatval_t E = 0.1;
    lbfgsfloatval_t Beta = 2.0;
    lbfgsfloatval_t EN[2];
    lbfgsfloatval_t* en;
    int TF = 0;
    
    en = &EN[0];
    int k = 0;
    int COUNT = 0;
    int m = 0;
    for (int i = 0; i < FIBER; i++) {
        m = sampling(E, Beta, en,i);
        
        F[i][0] = EN[0];
        F[i][1] = EN[1];
        if (m == 1) {
            break;
        }

    }

    for (int i = 0; i < fibersample.size(); i++) {
        fibersample[i][0] = F[i][0];
        fibersample[i][1] = F[i][1];
        fibersample[i][2] = 0.0;
        
        for (int h = 0; h < CIRCLE; h++) {
            lbfgsfloatval_t Gsin = sin((2.0 / CIRCLE) * h * pi);
            lbfgsfloatval_t Gcos = cos((2.0 / CIRCLE) * h * pi);
            fibersample2[h][i][0] = fibersample[i][0] + Gsin * Fradius;
            fibersample2[h][i][1] = fibersample[i][1] + Gcos * Fradius;
            fibersample2[h][i][2] = fibersample[i][2];
        }
        
    }
    for (int i = 0; i < Rayply; i++) {
        RayPly[i][0] = cos(2.0 * pi * i / Rayply);
        RayPly[i][1] = sin(2.0 * pi * i / Rayply);
        RayPly[i][2] = 0.0;
        //printf("x %f y %f z %f\n", RayPly[i][0], RayPly[i][1], RayPly[i][2]);
    }
     
    GLfloat A=0.0, B=0.0, C=0.0, D = 0.0;
    GLfloat t[2][FIBER][Rayply];
    for (int j = 0; j < Rayply; j++) {
        GLfloat tmp = 1.0;
        int tmpind = 360;
        int tmpcnt = 2;
        int ikeep = 0;
        
        for (int i = 0; i < fibersample.size(); i++) {

            //方向ベクトルと繊維の判別式を求める
            /*

            (x-a)^2 + (y-b)^2 = r^2
            (c - ct -a)^2 + (s - st - b)^2 = r^2
            (c^2 - 2cct - 2ac + cctt + 2act  + a^2) + (s^2 - 2sst - 2bs + sstt + 2bst  + b^2) = r^2

            (c^2 + s^2)*t^2  -  (2ac+2bs-2c^2-2s^2)*t  +  (c^2-2ac+a^2+s^2-2bs+b^2-r^2) = 0
              
            */
            /*
            A = RayPly[j][0] * RayPly[j][0] + RayPly[j][1] * RayPly[j][1];
            B = -(2 * fibersample[i][0] * RayPly[j][0] + 2 * fibersample[i][1] * RayPly[j][1] - 2 * RayPly[j][0] * RayPly[j][0] - 2 * RayPly[j][1] * RayPly[j][1]);
            C = RayPly[j][0] * RayPly[j][0] - 2 * fibersample[i][0] * RayPly[j][0] + fibersample[i][0] * fibersample[i][0]
                + RayPly[j][1] * RayPly[j][1] - 2 * fibersample[i][1] * RayPly[j][1] + fibersample[i][1] * fibersample[i][1] - Fradius * Fradius;
            */
            A = 1.0 + (RayPly[j][1] * RayPly[j][1]) / ((RayPly[j][0] + 1e-12 )* (RayPly[j][0] + 1e-12));
            B = -2.0 * (fibersample[i][0] + (RayPly[j][1] / (RayPly[j][0] + 1e-12)) * fibersample[i][1]);
            C = fibersample[i][0] * fibersample[i][0] + fibersample[i][1] * fibersample[i][1] - Fradius * Fradius;
            // 判別式をする A*t^2 + B*t + C = 0
            D = B * B - 4 * A * C;
            //printf("A = %f\nB = %f\nC = %f\nD = %f\n", A, B, C, D);
            t[0][i][j] = 1.0;
            t[1][i][j] = 1.0;
            Visible[j][0] = 0.0;
            Visible[j][1] = 0.0;
            Visible[j][2] = 0.0;
            if (D >= 0) {
                t[0][i][j] = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
                t[1][i][j] = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
                GLfloat tmp0 = t[0][i][j];// *RayPly[j][0];
                GLfloat tmp1 = (RayPly[j][1] / (RayPly[j][0] + 1e-12)) * t[0][i][j];// *RayPly[j][1];
                if (tmp > sqrt((RayPly[j][0] - tmp0) * (RayPly[j][0] - tmp0) + (RayPly[j][1] - tmp1) * (RayPly[j][1] - tmp1))) {
                    tmp = sqrt((RayPly[j][0] - tmp0) * (RayPly[j][0] - tmp0) + (RayPly[j][1] - tmp1) * (RayPly[j][1] - tmp1));
                    tmpind = i;
                    tmpcnt = 0;
                }
                GLfloat tmp2 = t[1][i][j];// *RayPly[j][0];
                GLfloat tmp3 = (RayPly[j][1] / (RayPly[j][0] + 1e-12)) * t[1][i][j];// *RayPly[j][1];
                if (tmp > sqrt((RayPly[j][0] - tmp2) * (RayPly[j][0] - tmp2) + (RayPly[j][1] - tmp3) * (RayPly[j][1] - tmp3))) {
                    tmp = sqrt((RayPly[j][0] - tmp2) * (RayPly[j][0] - tmp2) + (RayPly[j][1] - tmp3) * (RayPly[j][1] - tmp3));
                    tmpind = i;
                    tmpcnt = 1;
                }

                //printf("                   %f %f\n", t[0][i][j], t[1][i][j]);
                //printf("%f %f \n", Visible[0][i][j], Visible[1][i][j]);
            }

        }
        if (tmpind < 360 && tmp > 0.0 && tmpcnt < 2) {
            Visible[j][0] = t[tmpcnt][tmpind][j];
            Visible[j][1] = (RayPly[j][1] / (RayPly[j][0] + 1e-12)) * t[tmpcnt][tmpind][j];
            Visible[j][2] = 0.0;
            for (int k = 0; k < 3; k++) {
                PLYnormal[j][k] = Normal(fibersample, Visible, tmpind, j, k);
            }

        }
    }

    Shadowing();
    for (int j = 0; j < length; j++) {
        for (int i = 0; i < Rayply; i++) {
            float th = ((((float)j) * 2.0 * pi) / (float)plypitch);// +((float)i / (float)Rayply) * 2.0 * pi);//zはbθで、今回zはjのことだから、これを使ってθを求める
            RayPly10[j][i][0] = RayPly[i][0]; //Rradius *cos(th);
            RayPly10[j][i][1] = RayPly[i][1]; //Rradius *sin(th);
            RayPly10[j][i][2] = (j - 5.0) * pi;
            tex1D10[j][i][0] = (cos(th) * PLYnormal[i][0] - sin(th) * PLYnormal[i][1]);// *(Shadow[i] / 180.0);
            tex1D10[j][i][1] = (sin(th) * PLYnormal[i][0] + cos(th) * PLYnormal[i][1]);// *(Shadow[i] / 180.0);
            tex1D10[j][i][2] = PLYnormal[i][2];// *(Shadow[i] / 180.0);
            /*
            RayPly2[i][0] = RayPly[i][0];
            RayPly2[i][1] = RayPly[i][1];
            RayPly2[i][2] = -5.0;
            tex1D[i][0] = PLYnormal[i][0];
            tex1D[i][1] = PLYnormal[i][1];
            tex1D[i][2] = PLYnormal[i][2];
            */
            /*
            if (Shadow[i] < 60) {
                tex1D[i][0] = 0.0;
                tex1D[i][1] = 0.0;
                tex1D[i][2] = 0.0;
            }
            */
            
            
        }
        //printf("[%d] %f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n\n",j , tex1D10[j][0][0], tex1D10[j][0][1], tex1D10[j][0][2], tex1D10[j][89][0], tex1D10[j][89][1], tex1D10[j][89][2], tex1D10[j][179][0], tex1D10[j][179][1], tex1D10[j][179][2], tex1D10[j][269][0], tex1D10[j][269][1], tex1D10[j][269][2]);
    }
    
    
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 50);
    glutInitWindowSize(300, 300);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Hello GLUT!!");
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    init();
    glutMainLoop();
    return 0;
}
