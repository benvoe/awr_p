#include <iostream>
#include "Risiko_Spiel.h"

int main() {

    // AUSWAHL DER PRAKTIKUMSAUFGABEN

    bool p1_1  = false;
    bool p1_2a = false;
    bool p1_2b = false;
    bool p1_3a = false;
    bool p1_3b = false;
    bool p1_4  = false;

    bool p2_1a = false;
    bool p2_1b = false;
    bool p2_2a = false;
    bool p2_2b = false;

    bool p3 = false;

    bool p4 = false;

    bool p5 = false;

    bool p6 = true;

    CMyMath math = CMyMath();

    //cout << math.NueberK(32,0) << endl;

    //cout << math.berechneK0( 50, 0.2, 0.05) << endl;
    //cout << math.berechneK1( 50, 0.2, 0.05) << endl;

    //cout << math.integral_e(5) << endl;

    /*
    for(int i = 1; i <= 4; i++){
        //cout << "Integral " << i << ": " << 1/sqrt( 2 * M_PI ) * ( math.F(i) - math.F(-50) ) << endl;
    }
    */

    /*
    vector<double> res = vector<double>(10,0);
    int poisson = 0;
    int n = 100000;

    for( int i = 0; i < n; i++ ){
        poisson = math.generierePoisson( 2, 3.7);

        if(poisson < res.size()){
            res[poisson]++;
        }
    }

    for(int i = 0; i < res.size(); i++){
        res[i] /= n;
        cout << res[i] << endl;
    };
    */

    if(true) {

        cout << "Mode: 0" << endl;
        for (int i = 1; i < 10; i++) {
            cout << i / 10.0 << " : ";
            math.mittlereAufenthaltszeit(1000000, i / 10.0, 1, 0);
        }
        cout << endl;
        cout << "Mode: 1" << endl;

        for (int i = 1; i < 10; i++) {
            cout << i / 10.0 << " : ";
            math.mittlereAufenthaltszeit(1000000, i / 10.0, 1, 1);
        }
        cout << endl;
        cout << "Mode: 2" << endl;

        for (int i = 1; i < 10; i++) {
            cout << i / 10.0 << " : ";
            math.mittlereAufenthaltszeit(1000000, i / 10.0, 1, 2);
        }

        for (int i = 0; i < 10; i++) {
            //cout << math.natuerlicheZahlinN_0(2) << endl;
        }
    }


    if(p5) {

        // AUFGABE 5-1

        cout << "AUFGABE 1A:" << endl;
        math.berechneAkzeptanzbereich(50, 0.2, 0.9, 0);
        cout << endl;
        math.berechneAkzeptanzbereich(50, 0.2, 0.9, 1);
        cout << endl;
        math.berechneAkzeptanzbereich(50, 0.2, 0.9, -1);
        cout << endl << endl;


        // AUFGABE 5-2

        cout << "AUFGABE 2:" << endl;

        int counter = 0;
        int N = 100000;

        for( int i = 0; i < N; i++){
            counter += math.hypothesentest();
        }

        cout << "Hypothese liegt zu " << double(counter)/ double(N) << "% im Akzeptanzbereich." << endl;

    }

    // AUFGABE 6-1




    // PRAKTIKUM 1

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p1_1){
        math.wahrscheinlichkeitRisiko(3,2);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p1_2a){
        int n = 6;
        double ziehungen = 2000000.00;
        vector<int> ergebnisse = vector<int>(n, 0);

        for (int i = 0; i < ziehungen; i++) {
            int x = math.natuerlicheZahlinN(n);
            ergebnisse[x-1]++;
        }

        for(int i = 0; i < n; i++){
            std::cout <<  i+1 << " gezogen " << ergebnisse[i] << " mal" << " ( " << ergebnisse[i] / ziehungen << " % ) " << std::endl;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p1_2b){
        int in_p = 0, nicht_p = 0;
        double p = 0.25;
        double ziehungen = 2000000;

        for (int i = 0; i < ziehungen; i++) {
            if (math.bernoulliSimulation(p)) {
                in_p++;
            } else
                nicht_p++;
        }

        std::cout << "p gezogen " << in_p << " mal" << " ( " << in_p / ziehungen << " ) " << std::endl;
        std::cout << "1-p gezogen " << nicht_p << " mal" << " ( " << nicht_p / ziehungen << " ) " << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p1_3a){
        Risiko_Spiel spiel = Risiko_Spiel();

        spiel.spielen(3, 2, 2000000);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p1_3b){
        math.streuungBitfehler(0.01, 2, 20000, 2000);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p1_4){
        if(false){
            for(int i = 1; i <= 15; i += 1){
                math.streuungRisiko(3, 2, i, 100000, false);
            }
        }
        else {
            for(int i = 1; i <= 50; i += 1){
                math.streuungBitfehler(0.01, 2, i, 1000, false);
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // PRAKTIKUM 2

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p2_1a){ //  int iterationen, int a, int b, int anz, double lambda
        math.generiereHistogrammExponentialverteilt( 5000000, 0, 5, 50, 2 );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p2_1b){ //  int iterationen, int dx, int anz, int mittel, double varianz
        math.generiereHistogrammNormalverteilt( 5000000, 3, 50, 0, 0.5 );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p2_2a){ //  int b, double lambda, int bereich, double genauigkeit, int iterationen
        math.simuliereOptimaleWartezeitExponentiell( 10, 1, 5, 0.1, 100000);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(p2_2b){ //  int b, double stdAbweichung, int bereich, double genauigkeit, int iterationen
        math.simuliereOptimaleWartezeitNormalverteilt( 10, 1, 2, 0.1, 1000000 );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // PRAKTIKUM 3

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( p3 ){

        //Normalverteilt
        cout << "Mittelwerte - Normalverteilung" << endl;
        math.stichprobemMittelNormalverteilung(5000000, 10, 0, 1);
        math.stichprobemMittelNormalverteilung(5000000, 10, 0, 1.5);
        math.stichprobemMittelNormalverteilung(5000000, 10, 0, 2);
        math.stichprobemMittelNormalverteilung(5000000, 10, 0, 2.5);
        math.stichprobemMittelNormalverteilung(5000000, 10, 0, 3);

        // Exponentialverteilt
        cout << "Mittelwerte - Exponentialverteilung" << endl;
        math.stichprobemMittelExponentialverteilung(5000000, 10, 0.5);
        math.stichprobemMittelExponentialverteilung(5000000, 10, 1);
        math.stichprobemMittelExponentialverteilung(5000000, 10, 1.5);
        math.stichprobemMittelExponentialverteilung(5000000, 10, 2);
        math.stichprobemMittelExponentialverteilung(5000000, 10, 2.5);
        math.stichprobemMittelExponentialverteilung(5000000, 10, 3);


        cout << "Mittelwerte - Exponentialverteilung (variables n)" << endl;
        math.stichprobemMittelExponentialverteilung(5000, 100000, 0.5);
        math.stichprobemMittelExponentialverteilung(5000, 10000, 0.5);
        math.stichprobemMittelExponentialverteilung(5000, 1000, 0.5);
        math.stichprobemMittelExponentialverteilung(5000, 100, 0.5);
        math.stichprobemMittelExponentialverteilung(5000, 10, 0.5);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // PRAKTIKUM 4

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( p4 ){

        int n = 1000;

        math.simuliereKonfidenzbereich( 200000, n );

        for( int N = 10; N <= 100000; N = N * 10) {
            //n = natuerlicheZahlinN( 100 ); cout << n << endl;
            cout << "\nAnzahl N = " << N << endl;
            for(int i = 0; i < 5; i++ ) {
                cout << "Anzahl n = " << n << " - " ;
                math.simuliereKonfidenzbereich(N, n);
            }
        }

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return 0;
}