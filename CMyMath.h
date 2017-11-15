//
// Created by Benedikt Völker on 25.04.17.
//

#ifndef AWR_PRAKTIKUM_CMYMATH_H
#define AWR_PRAKTIKUM_CMYMATH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
using namespace std;



struct mittelwert_varianz {
    double mittelwert = 0;
    double varianz = 0;
};

struct konfidenzbereich {
    double begin = 0;
    double end = 0;
};

struct akzeptanzbereich {
    double k0 = 0;
    double k1 = 0;
};

struct kunde {
    double t0 = 0;
    double t1 = 0;
};

class CMyMath {
public:
    CMyMath();

    // Hilfsfunktionen
    string cs( double zahl );
    double integral_e ( double x );
    double F( double x );

    double NueberK( int n, int k );
    int fakultaet( int n );

    int natuerlicheZahlinN ( int n );
    int natuerlicheZahlinN_0 ( int n );
    int natuerlicheZahlinN_neu ( int n );
    int natuerlicheZahlinAB( int a, int b );

    bool bernoulliSimulation ( double p );

    double zufallszahl_in_0_1 ();
    double generiereExponentialverteilung ( double lambda );
    double generiereNormalverteilung ( double mittelwert, double varianz );

    // PRAKTIKUM 1
    void wahrscheinlichkeitRisiko(int k, int l);

    double wahrscheinlichkeitBitfehlerN32( double p, int f, int n );
    void streuungBitfehler( double p, int f, int n, int N, bool ausgabe = true);
    void streuungRisiko( int k, int l, int n, int N, bool ausgabe = true);

    // PRAKTIKUM 2
    void generiereHistogrammExponentialverteilt ( int iterationen, int a, int b, int anz, double lambda );
    void generiereHistogrammNormalverteilt ( int iterationen, int dx, int anz, int mittel, double varianz );

    // Exponentialverteilung
    void simuliereOptimaleWartezeitExponentiell ( int b, double lambda, int bereich, double genauigkeit, int iterationen );
    double wartezeitExponentiell( int b, double x, double lambda );

    // Normalverteilung
    void simuliereOptimaleWartezeitNormalverteilt ( int b, double stdAbweichung, int bereich, double genauigkeit, int iterationen );
    double wartezeitNormalverteilt( int b, double x, double stdAbweichung );

    // PRAKTIKUM 3
    mittelwert_varianz stichprobeNormalverteilung( int n, double erwartungswert, double standardabweichug );
    void stichprobemMittelNormalverteilung(int N, int n, double erwartungswert, double standardabweichung );

    double stichprobeExponentialverteilung ( int n, double lambda );
    void stichprobemMittelExponentialverteilung (int N, int n, double lambda );

    // PRAKTIKUM 4
    void simuliereKonfidenzbereich( int N, int n);
    bool stichprobeKonfidenzbereich( int n, double p );
    konfidenzbereich berechneKonfidenzbereich( int n, int k );

    // PRAKTIKUM 5
    akzeptanzbereich berechneAkzeptanzbereich( int n, double p, double a, int mode);
    double berechneK0( int n, double p, double a );
    double berechneK1( int n, double p, double a );

    bool hypothesentest();

    // PRAKTIKUM 6
    int generierePoisson( double a, double t = 1 );
    void mittlereAufenthaltszeit( int N, double a, double b, int mode );

};


#endif //AWR_PRAKTIKUM_CMYMATH_H
