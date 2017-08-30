//
// Copyright(c) 2017 by Nobuo NAKAGAWA @ Polyphony Digital Inc.
//
// We're Hiring!
// http://www.polyphony.co.jp/recruit/
//
#include <cstdlib>
#include <cstdint>
#include <vector>
#include "glm/glm.hpp"
#if defined(WIN32)
#include <GL/glut.h>
#ifndef _DEBUG
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
#endif // _DEBUG
#elif defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#endif // MACOSX

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct sDistanceConstraint {
  int32_t p1;
  int32_t p2;
  GLfloat rest;
  GLfloat k;
  GLfloat k_prime;
};

struct sTriangleConstraint {
  int32_t p1;
  int32_t p2;
  int32_t p3;
};

struct sBendingConstraint {
  int32_t p1;
  int32_t p2;
  int32_t p3;
  int32_t p4;
  GLfloat k_prime;
};

struct sParticle {
 glm::vec3    P; // position
 glm::vec3   pP; // predicted position
 glm::vec3   dP; // position correction
 glm::vec2   UV; // texture cordinate
 glm::vec3    V; // velocity
 glm::vec3    F; // force
 GLfloat   invW; // inverse_mass(=1.0f / mass)
 glm::vec3   Ri; // Ri = Xi - Xcm
};

struct sParam {
  float Sxx;
  float Syy;
  float Sxy;
  sParam() : Sxx(1.0f),
             Syy(1.0f),
             Sxy(1.0f) {
  }
};

struct sApplication {
  std::vector<sDistanceConstraint> d_constraints;
  std::vector<sTriangleConstraint> t_constraints;
  std::vector<sBendingConstraint>  b_constraints;
  std::vector<GLfloat>             phi0;          // initiali dihedral angle between adjacent triangles
  std::vector<GLushort>            indices;
  std::vector<glm::mat2>           Qinv;          // inverse of material coordinates
  std::vector<sParticle>           particles;
  sParam                           param;
};

static const glm::ivec2 CLOTH_QUADS(20, 20);
static const int32_t    TOTAL_POINTS = (CLOTH_QUADS.x + 1) * (CLOTH_QUADS.y + 1);
static const int32_t    SIZE  = 4;
static const GLfloat    HSIZE = SIZE / 2.0f;
static const GLfloat    MASS  = 1.0f;

static const GLfloat    STRETCH        = 0.5f;      // 0.0f-1.0f
static const GLfloat    BEND           = 0.75f;     // 0.0f-1.0f
static const GLfloat    DAMP           = 0.000125f; // 0.0f-1.0f
static const GLfloat    GLOBAL_DAMPING = 0.98f;     // 0.0f-1.0f
static const uint32_t   SOLVER_ITERATIONS = 10;
static const glm::vec3  GRAVITY           = glm::vec3(0.0f, -9.81f, 0.0f);

sApplication g_App;

void add_bending_constraint(int32_t pa, int32_t pb, int32_t pc, int32_t pd, float k) {
  sBendingConstraint c;
  c.p1 = pa;
  c.p2 = pb;
  c.p3 = pc;
  c.p4 = pd;
  c.k_prime = 1.0f - pow((1.0f - k), 1.0f / SOLVER_ITERATIONS);
  if (c.k_prime > 1.0)
    c.k_prime = 1.0f;
  g_App.b_constraints.push_back(c);
}

glm::vec3 get_normal(int ind0, int ind1, int ind2) {
  std::vector<sParticle>& p = g_App.particles;
  glm::vec3 e1 = p[ind0].P - p[ind1].P;
  glm::vec3 e2 = p[ind2].P - p[ind1].P;
  return glm::normalize(glm::cross(e1,e2));
}

float get_dihedral_angle(sBendingConstraint c, float& d, glm::vec3& n1, glm::vec3& n2) {
  n1 = get_normal(c.p1, c.p2, c.p3);
  n2 = get_normal(c.p1, c.p2, c.p4); 
  d = glm::dot(n1, n2);
  return acos(d);
} 

void add_triangle_constraint(int pa, int pb, int pc, float k) {
  sTriangleConstraint c;
  c.p1 = pa;
  c.p2 = pb;
  c.p3 = pc;
  g_App.t_constraints.push_back(c);
}

void init(int argc, char* argv[]) {
  glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
  glEnable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  g_App.indices.resize(CLOTH_QUADS.x * CLOTH_QUADS.y * 2 * 3);
  g_App.Qinv.resize(g_App.indices.size() / 3);
  g_App.particles.resize(TOTAL_POINTS);

  int cnt = 0;
  std::vector<sParticle>& p = g_App.particles;
  for(int32_t j = 0; j <= CLOTH_QUADS.y; j++) { // fill in positions
    for(int32_t i = 0; i <= CLOTH_QUADS.x; i++) {
      float u = (float(i) / CLOTH_QUADS.x);
      float v = (float(j) / CLOTH_QUADS.y);
      p[cnt].UV = glm::vec2(u, v);
      p[cnt].P  = p[cnt].pP = glm::vec3((u * 2 - 1) * HSIZE, SIZE + 1.0f, v * SIZE);
      p[cnt].invW  = 1.0f / MASS;
      cnt++;
    }
  }
  p[0].invW = p[CLOTH_QUADS.x - 1].invW = 0.0f;      // fixed points
  GLushort* id = &(g_App.indices[0]);
  for(int32_t i = 0; i < CLOTH_QUADS.y; i++) { // fill in indices
    for(int32_t j = 0; j < CLOTH_QUADS.x; j++) {
      int i0 = i  * (CLOTH_QUADS.x+1) + j;     // i0 - i1
      int i1 = i0 + 1;                         //  |   |
      int i2 = i0 + (CLOTH_QUADS.x+1);         // i2 - i3
      int i3 = i2 + 1;
      if ((j + i) % 2) {                       // index
        *id++ = i0; *id++ = i2; *id++ = i1;    // 0-2-1
        *id++ = i1; *id++ = i2; *id++ = i3;    // 1-2-3
      } else {
        *id++ = i0; *id++ = i2; *id++ = i3;    // 0-2-3
        *id++ = i0; *id++ = i3; *id++ = i1;    // 0-3-1
      }
    }
  }
  std::vector<GLushort>& indices = g_App.indices;
  int32_t indices_size = indices.size();
  int32_t count        = 0;
  for(int32_t i = 0; i < indices_size; i+= 3) { // Material Coordinates for Triangles Eq.37(p.9)
    int i0 = indices[i  ];
    int i1 = indices[i+1];
    int i2 = indices[i+2];
    glm::vec3   x0   = p[i0].P;
    glm::vec3   x1   = p[i1].P;
    glm::vec3   x2   = p[i2].P;
    glm::vec2   u0   = p[i0].UV;
    glm::vec2   u1   = p[i1].UV;
    glm::vec2   u2   = p[i2].UV;
    glm::mat2x3 P    = glm::mat2x3((x1-x0), (x2-x0));
    glm::mat2x2 U    = glm::mat2x2((u1-u0), (u2-u0));
    glm::mat2x2 Uinv = glm::inverse(U);
    glm::mat2x3 T    = P * Uinv;
    glm::vec3   n1(T[0]);
    glm::vec3   n2(T[1]);
    n1 = glm::normalize(n1);
    n2 = glm::normalize(n2);
    glm::mat2x2 C     = glm::transpose(glm::mat2x3(n1, n2)) * P;
    g_App.Qinv[count] = glm::inverse(C);
    count++;
  }
  // setup bending constraints
  int v = CLOTH_QUADS.y + 1;
  int u = CLOTH_QUADS.x + 1;
  for(int32_t i = 0; i < CLOTH_QUADS.y; i++) {
    for(int32_t j = 0; j < CLOTH_QUADS.x; j++) {
      int p1 = i * (CLOTH_QUADS.x + 1) + j;
      int p2 = p1 + 1;
      int p3 = p1 + CLOTH_QUADS.x + 1;
      int p4 = p3 + 1;
      if ((j + 1) % 2) {
        add_bending_constraint(p3, p2, p1, p4, BEND);
      } else {
        add_bending_constraint(p4, p1, p3, p2, BEND);
      }
    }
  }
  GLfloat d;
  glm::vec3 n1, n2;
  int32_t b_constraints_size = g_App.b_constraints.size();
  g_App.phi0.resize(b_constraints_size);
  for(int32_t i = 0 ; i < b_constraints_size; i++) {
    g_App.phi0[i] = get_dihedral_angle(g_App.b_constraints[i], d, n1, n2);
  }
  for(int32_t i = 0; i < indices_size; i+=3) { // triangle constraints
    int32_t i0 = indices[i    ];
    int32_t i1 = indices[i + 1];
    int32_t i2 = indices[i + 2];
    add_triangle_constraint(i0, i1, i2, STRETCH);
  }
  
}

void display(void){
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(  0.0f,  // pos
             -5.0f,
             -8.0f,
              0.0f,  // tgt
              0.0f,
              0.0f,
              0.0f,
              1.0f,
              0.0f); // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_NORMALIZE);
#if 0
  glFrontFace(GL_CW);
  glutSolidTeapot(1.0f);
  glFrontFace(GL_CCW);
#else
  if (1) {
    glFrontFace(GL_CW);
  }
  glBegin(GL_TRIANGLES);
  glColor3f(1.0f, 1.0f, 1.0f);
  std::vector<GLushort>&  indices      = g_App.indices;
  int32_t                 indices_size = indices.size();
  std::vector<sParticle>& particles    = g_App.particles;
  for(int32_t i = 0; i < indices_size; i+=3) {
    int32_t   i0  = indices[i  ];
    int32_t   i1  = indices[i+1];
    int32_t   i2  = indices[i+2];
    glm::vec3 p1  = particles[i0].P;
    glm::vec3 p2  = particles[i1].P;
    glm::vec3 p3  = particles[i2].P;
    glm::vec2 uv1 = particles[i0].UV;
    glm::vec2 uv2 = particles[i1].UV;
    glm::vec2 uv3 = particles[i2].UV;
    glm::vec3 N = glm::normalize(glm::cross((p2-p1), (p3-p1)));
    glNormal3f(N.x, N.y, N.z);
    glVertex3f(p1.x, p1.y, p1.z);
    glVertex3f(p2.x, p2.y, p2.z);
    glVertex3f(p3.x, p3.y, p3.z);
  }
  glEnd();

  if (1) {
    glDisable(GL_LIGHTING);
    glPointSize(5);
    glBegin(GL_POINTS);
    for(int32_t i = 0; i < TOTAL_POINTS; i++) {
      glm::vec3 p = particles[i].P;
      glColor3f(1.0f, 0.0f, 0.0f);
      glVertex3f(p.x, p.y, p.z);
    }
    glEnd();
    glEnable(GL_LIGHTING);
  }
#endif
  glutSwapBuffers();
}

void reshape(int width, int height){
  static GLfloat lightPosition[4] = {0.0f, 250.0f, 55.0f, 1.0f};
  static GLfloat lightDiffuse[3]  = {0.5f,  0.5f,  0.5f      };
  static GLfloat lightAmbient[3]  = {0.25f, 0.25f, 0.25f     };
  static GLfloat lightSpecular[3] = {1.0f,  1.0f,  1.0f      };

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glViewport(0, 0, width, height);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  GLdouble fovy = 45.0;
  gluPerspective(fovy, (double)width / (double)height, 1.0, 10000.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDiffuse);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  lightAmbient);
  glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
}

glm::mat3x2 outer_product(const glm::vec3& a, const glm::vec2& b) {
  glm::mat3x2 ret;
  ret[0][0] = a.x * b.x;
  ret[1][0] = a.y * b.x;
  ret[2][0] = a.z * b.x;
  ret[0][1] = a.x * b.y;
  ret[1][1] = a.y * b.y;
  ret[2][1] = a.z * b.y;

  return ret;
}

void compute_forces() {
  std::vector<sParticle>& particles = g_App.particles;
  for(int32_t i = 0; i < TOTAL_POINTS; i++) {
    if (particles[i].invW > 0.0f) {     // invW = 0.0f means fixed points
      particles[i].F  = glm::vec3(0);
      particles[i].F += GRAVITY * MASS; // F = ma
    }
  }
}

void integrate_explicit(GLfloat dt) {
  std::vector<sParticle>& particles = g_App.particles;
  for(int32_t i = 0; i < TOTAL_POINTS; i++) {
#if 1
    particles[i].V *= GLOBAL_DAMPING;
    particles[i].V += particles[i].F * dt * particles[i].invW;
#endif
    if (particles[i].invW < 0.0f) {     // invW = 0.0f means fixed points
      particles[i].pP  = particles[i].P; // fixed
    } else {
      particles[i].pP += particles[i].V * dt;
    }
  }
}

void update_triangle_constraint() {
  std::vector<sTriangleConstraint>& t_constraints      = g_App.t_constraints;
  std::vector<sParticle>&           particles          = g_App.particles;
  int32_t                           t_constraints_size = t_constraints.size();
  for(int32_t i = 0; i < t_constraints_size; i++) {
    sTriangleConstraint c = t_constraints[i];
    float w1 = particles[c.p1].invW;
    float w2 = particles[c.p2].invW;
    float w3 = particles[c.p3].invW;
    glm::vec3 x0 = particles[c.p1].pP;
    glm::vec3 x1 = particles[c.p2].pP;
    glm::vec3 x2 = particles[c.p3].pP;
    glm::mat2x2 C = g_App.Qinv[i];
    glm::vec2 C1 = C[0];
    glm::vec2 C2 = C[1];
    glm::mat2x3 P =(glm::mat2x3((x1-x0),(x2-x0)));
    glm::mat2x3 F = P * C;
    glm::vec3 f1(F[0]);
    glm::vec3 f2(F[1]);
    float S11 = glm::dot(f1, f1);                              // Eq.28
    float S12 = glm::dot(f1, f2);
    float S22 = glm::dot(f2, f2);
    glm::mat3x2 mat11 = (2.0f * outer_product(f1, C1));        // Eq.29
    glm::mat3x2 mat12 = (outer_product(f2, C1) + outer_product(f1, C2));
    glm::mat3x2 mat22 = (2.0f * outer_product(f2, C2));
    glm::vec3 ds1_dp1(mat11[0][0], mat11[1][0], mat11[2][0]);  // stretch constraints
    glm::vec3 ds1_dp2(mat11[0][1], mat11[1][1], mat11[2][1]);
    glm::vec3 ds1_dp0 = -(ds1_dp1 + ds1_dp2);
    glm::vec3 ds2_dp1(mat22[0][0], mat22[1][0], mat22[2][0]);
    glm::vec3 ds2_dp2(mat22[0][1], mat22[1][1], mat22[2][1]);
    glm::vec3 ds2_dp0 = -(ds2_dp1 + ds2_dp2);
    glm::vec3 ds12_dp1(mat12[0][0], mat12[1][0], mat12[2][0]); // shear constraints
    glm::vec3 ds12_dp2(mat12[0][1], mat12[1][1], mat12[2][1]);
    glm::vec3 ds12_dp0 = -(ds12_dp1 + ds12_dp2);
    glm::vec3 dp0(0,0,0);
    glm::vec3 dp1(0,0,0);
    glm::vec3 dp2(0,0,0);
    float weight_sum11 = w1 * glm::dot(ds1_dp0, ds1_dp0) +     // 1) stretch - S11(Eq.33)
                         w2 * glm::dot(ds1_dp1, ds1_dp1) +
                         w3 * glm::dot(ds1_dp2, ds1_dp2);
    float lambda11 = (S11 - 1) / weight_sum11;
    S11 = -g_App.param.Sxx * lambda11;
    dp0 += w1 * S11 * ds1_dp0;
    dp1 += w2 * S11 * ds1_dp1;
    dp2 += w3 * S11 * ds1_dp2;
    float weight_sum22 = w1 * glm::dot(ds2_dp0, ds2_dp0) +     // 2) stretch - S22(Eq.33)
                         w2 * glm::dot(ds2_dp1, ds2_dp1) +
                         w3 * glm::dot(ds2_dp2, ds2_dp2);
    float lambda22 = (S22 - 1) / weight_sum22;
    S22 = -g_App.param.Syy * lambda22;
    dp0 += w1 * S22 * ds2_dp0;
    dp1 += w2 * S22 * ds2_dp1;
    dp2 += w3 * S22 * ds2_dp2;
    float weight_sum12 = w1 * glm::dot(ds12_dp0, ds12_dp0) +   // 4) Shear - S12(Eq.34)
                         w2 * glm::dot(ds12_dp1, ds12_dp1) +
                         w3 * glm::dot(ds12_dp2, ds12_dp2);
    float lambda12 = S12 / weight_sum12;
    S12 = -g_App.param.Sxy * lambda12;
    dp0 += w1 * S12 * ds12_dp0;
    dp1 += w2 * S12 * ds12_dp1;
    dp2 += w3 * S12 * ds12_dp2;
    if (w1 > FLT_EPSILON) {
      particles[c.p1].dP += dp0;
      particles[c.p1].pP += dp0;
    }
    if (w2 > FLT_EPSILON) {
      particles[c.p2].dP += dp1;
      particles[c.p2].pP += dp1;
    }
    if (w3 > FLT_EPSILON) {
      particles[c.p3].dP += dp2;
      particles[c.p3].pP += dp2;
    }
  }
}

void update_bending_constraint() {
  std::vector<sBendingConstraint>& b_constraints      = g_App.b_constraints;
  std::vector<sParticle>&          particles          = g_App.particles;
  int32_t                          b_constraints_size = b_constraints.size();
  for(int32_t i = 0; i < b_constraints_size; i++) {
    sBendingConstraint c = b_constraints[i];
    glm::vec3 p1 = particles[c.p1].pP;
    glm::vec3 p2 = particles[c.p2].pP - p1;
    glm::vec3 p3 = particles[c.p3].pP - p1;
    glm::vec3 p4 = particles[c.p4].pP - p1;
    glm::vec3 p3p1 = p3 - p1;
    glm::vec3 p4p1 = p4 - p1;
    glm::vec3 p3p2 = p3 - p2;
    glm::vec3 p4p2 = p4 - p2;
    glm::vec3 e  = p4 - p3;          // prepare e, n1 and n2 using Eq.44, 45, 46
    glm::vec3 n1 = glm::cross(p3p1, p4p1);
    float denN1  = glm::dot(n1, n1); // glm::length(n1) * glm::length(n1);
    if (denN1 <= FLT_EPSILON)
      return;
    n1 = n1 / denN1;
    glm::vec3 n2 = glm::cross(p4p2, p3p2);
    float  denN2 = glm::dot(n2, n2);
    if (denN2 <= FLT_EPSILON)
      return;
    n2 = n2 / denN2;
    float d    = glm::dot(n1, n2);
    float phi  = acos(d);
    float sign = (glm::dot(glm::cross(n1, n2), e) > 0.0f) ? -1.0f : 1.0f;
    if      (d < -1.0f) { d = -1.0f; } // d = clamp(d, -1.0f, +1.0f);
    else if (d > +1.0f) { d = +1.0f; }
    if (d == -1.0f) { // acos(-1.0f) == M_PI
      phi = M_PI;
      if (phi == g_App.phi0[i])
        return;
      return; // TODO
    }
    if (d == 1.0f) {  // 180(triangles are planar)
      phi = 0.0f;
      if (phi == g_App.phi0[i])
        return;
    }
    float i_d = sqrt(1.0f - (d * d)) * (phi - g_App.phi0[i]);
    float lenE = glm::length(e);
    if (lenE <= FLT_EPSILON)
      return;
    glm::vec3 p1p4 = -p4p1;
    glm::vec3 p2p4 = -p4p2;
    glm::vec3 q1 = lenE * n1 * sign; // Eq. 40
    glm::vec3 q2 = lenE * n2 * sign; // Eq. 41
    glm::vec3 q3 = ((glm::dot(p1p4, e) / lenE) * n1 + (glm::dot(p2p4, e) / lenE) * n2) * sign; // Eq. 42
    glm::vec3 q4 = ((glm::dot(p3p1, e) / lenE) * n1 + (glm::dot(p3p2, e) / lenE) * n2) * sign; // Eq. 43
    float q1_len2 = glm::dot(q1, q1);
    float q2_len2 = glm::dot(q2, q2);
    float q3_len2 = glm::dot(q3, q3);
    float q4_len2 = glm::dot(q4, q4);
    float sum = particles[c.p1].invW * (q1_len2) + particles[c.p2].invW * (q2_len2) + particles[c.p3].invW * (q3_len2) + particles[c.p4].invW * (q4_len2);
    if (sum <= FLT_EPSILON)
      return;
    glm::vec3 dP1 = -( (particles[c.p1].invW * i_d) / sum) * q1;
    glm::vec3 dP2 = -( (particles[c.p2].invW * i_d) / sum) * q2;
    glm::vec3 dP3 = -( (particles[c.p3].invW * i_d) / sum) * q3;
    glm::vec3 dP4 = -( (particles[c.p4].invW * i_d) / sum) * q4;
    if (particles[c.p1].invW > 0.0f) {
      particles[c.p1].pP += dP1 * c.k_prime;
    }
    if (particles[c.p2].invW > 0.0f) {
      particles[c.p2].pP += dP2 * c.k_prime;
    }
    if (particles[c.p3].invW > 0.0f) {
      particles[c.p3].pP += dP3 * c.k_prime;
    }
    if (particles[c.p4].invW > 0.0f) {
      particles[c.p4].pP += dP4 * c.k_prime;
    }
  }
}

void integrate(GLfloat dt) {
  float inv_dt = 1.0f / dt;
  std::vector<sParticle>& particles = g_App.particles;
  for(int32_t i = 0; i < TOTAL_POINTS; i++) {
    particles[i].V = (particles[i].pP - particles[i].P) * inv_dt;
    particles[i].P = particles[i].pP;
  }
}

void idle(void) {
  compute_forces();
  integrate_explicit(1.0f / 60.0f);
  update_triangle_constraint();
  update_bending_constraint();
  integrate(1.0f / 60.0f);
  glutPostRedisplay();
}

void keyboard(unsigned char key , int x , int y) {
  switch(key){
  case 27:  exit(EXIT_SUCCESS); break; // ESC to quit
  }
}

int main(int argc, char* argv[]){
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowSize(640, 480);
  glutCreateWindow("Strain based dynamics");

  init(argc, argv);

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(idle);

  glutKeyboardFunc(keyboard);
/*
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutSpecialFunc(special);
*/
  glutMainLoop();
  return 0;
}
