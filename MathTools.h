//
// Created by Benedikt Völker on 11/04/2017.
//

#ifndef AWR_PRAKTIKUM_MATHTOOLS_H
#define AWR_PRAKTIKUM_MATHTOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

inline string cs ( double zahl ){
    string res = to_string(zahl);

    for(char& c : res) {
        if ( c == '.')
            c = ',';
    }

    return res;
}

inline int natuerlicheZahlinN ( int n ) {
    return rand() % n + 1;
}

inline int natuerlicheZahlinN_neu ( int n ) {
    int range = RAND_MAX / n;
    int i = 1;
    int random = rand();
    while( i * range < random ){
        i++;
        if(i > n){
            i = 1;
            random = rand();
        }
    }
    return i;
}

inline double zufallszahl_in_0_1 (){
    return double(rand()) / double(RAND_MAX);
}

inline bool bernoulliSimulation ( double p ) {
    return static_cast <double> (rand()) / static_cast <double> (RAND_MAX) <= p;
}

inline double generiereExponentialverteilung ( double lambda ){
    return - log( zufallszahl_in_0_1() ) / lambda;
}

inline double generiereNormalverteilung ( double mittelwert, double varianz ){
    return mittelwert + sqrt( -2 * varianz * log( zufallszahl_in_0_1() ) ) * sin( 2 * M_PI * zufallszahl_in_0_1() );
}

inline void generiereHistogrammExponentialverteilt ( int iterationen, int a, int b, int anz, double lambda ){
    double step = (b-a)/double(anz);
    double table[anz][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = a + (i+1) * step;
        table[i][1] = 0;
    }

    double zahl;
    for(int i = 0; i < iterationen; i++){
        zahl = generiereExponentialverteilung( lambda );
        for(int j = 0; j < anz; j++){
            if(table[j][0] > zahl ){
                table[j][1] += 1;
                break;
            }
        }
    }

    double scale = 1/table[0][1] * lambda;

    for(int i = 0; i < anz; i++){
        // cout << "Zwischen " << table[i][0] - step << " und " << table[i][0] << " : " << (table[i][1] * scale) << endl;
        cout << cs(table[i][1] * scale) << endl;
    }

}

inline void generiereHistogrammNormalverteilt ( int iterationen, int dx, int anz, int mittel, double varianz){
    double step = (2*dx)/double(anz);
    double table[50][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = - dx + (i+1) * step;
        table[i][1] = 0;
    }

    double zahl;
    for(int i = 0; i < iterationen; i++){
        zahl = generiereNormalverteilung(mittel,varianz);
        for(int j = 0; j < anz; j++){
            if(zahl > - dx && table[j][0] > zahl ){
                table[j][1] += 1;
                break;
            }
        }
    }

    int point = (int)(50 * ( dx+mittel ) / ( 2*dx ));

    double scale = 1/table[(int)(point)][1] / (sqrt(2*M_PI*varianz));

    for(int i = 0; i < anz; i++){
        // cout << "Zwischen " << table[i][0] - step << " und " << table[i][0] << " : " << (table[i][1] * scale) << endl;
        cout << cs(table[i][1] * scale) << endl;
    }
}

inline void simuliereOptimaleWartezeitExponentiell (int b, double lambda, int bereich, double genauigkeit ){
    int anz = bereich / genauigkeit;
    double table[anz][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = - bereich + (i+1) * genauigkeit;
        table[i][1] = 0;
    }

    double zahl;
    for(int i = 0; i < anz; i++){
        /*
        zahl = generiereNormalverteilung(mittel,varianz);
        for(int j = 0; j < anz; j++){
            if(zahl > - dx && table[j][0] > zahl ){
                table[j][1] += 1;
                break;
            }
        }
         */
    }

}

inline void wahrscheinlichkeitK3L2(){

    int zweifachSieg = 0, einfachSieg = 0, keinSieg = 0, moeglichkeiten = 0;
    vector<int> angreifer = {1, 1, 1};
    vector<int> verteidiger = {1, 1};
    vector<int> ang_sort;
    vector<int> ver_sort;

    while(angreifer.at(0) < 7) {
        while (verteidiger.at(0) < 7) {
            moeglichkeiten++;

            ang_sort = angreifer;
            sort(ang_sort.begin(), ang_sort.end(), std::greater<int>());

            ver_sort = verteidiger;
            sort(ver_sort.begin(), ver_sort.end(), std::greater<int>());


            if(ang_sort.at(0) > ver_sort.at(0) && ang_sort.at(1) > ver_sort.at(1)){
                zweifachSieg++;
            }
            else if(ang_sort.at(0) <= ver_sort.at(0) && ang_sort.at(1) <= ver_sort.at(1)){
                keinSieg++;
            }
            else {
                einfachSieg++;
            }

            verteidiger.at(1)++;
            if(verteidiger.at(1) > 6){
                verteidiger.at(1) = 1;
                verteidiger.at(0)++;
            }
            if(verteidiger.at(0) > 6){
                verteidiger.at(0) = 1;
                verteidiger.at(1) = 1;
                break;
            }
        }
        angreifer.at(2)++;
        if(angreifer.at(2) > 6){
            angreifer.at(2) = 1;
            angreifer.at(1)++;
        }
        if(angreifer.at(1) > 6){
            angreifer.at(1) = 1;
            angreifer.at(0)++;
        }
        if(angreifer.at(0) > 6){
            break;
        }
    }

    cout << "Wahrscheinlichkeit für doppelten Sieg: " << zweifachSieg/double(moeglichkeiten) << endl;
    cout << "Wahrscheinlichkeit für einfachen Sieg: " << einfachSieg/double(moeglichkeiten) << endl;
    cout << "Wahrscheinlichkeit für keinen Sieg: " << keinSieg/double(moeglichkeiten) << endl;
}

inline double wahrscheinlichkeitBitfehlerGroesser2UebertragungN32( double p, int n ){
    int mehrAls2 = 0;
    int bitfehler = 0;
    int wahrscheinlichkeit = 1;
    double result = 0;

    for(int i = 0; i < n; i++){
        for(int b = 0; b < 32; b++){
            if(bernoulliSimulation(p))
                bitfehler++;
        }
        if(bitfehler > 2){
            mehrAls2++;
        }
        bitfehler = 0;
    }
    result = mehrAls2 / double(n);

    // cout << "Die Wahrscheinlichkeit für mehr als 2 Bitfehler ist: " << result * 100 <<  "%" << endl;
    return result;
}

inline void streuungMehrAls2Bitfehler( double p, int n , int N){
    vector<float> stichproben = vector<float>();
    double mittelwert = 0, varianz = 0, standardAbweichung = 0;

    for(int i = 0; i < N; i++){
        stichproben.push_back( wahrscheinlichkeitBitfehlerGroesser2UebertragungN32( p, n ) );
        mittelwert += stichproben.back();
    }

    mittelwert /= n;

    for(int i = 0; i < n; i++){
        varianz += (stichproben.at(i) - mittelwert)*(stichproben.at(i) - mittelwert);
    }

    varianz /= n-1;

    standardAbweichung = sqrt(varianz);

    cout << "Mittelwert: " << mittelwert << endl;
    cout << "Varianz: " << varianz << endl;
    cout << "Die Standartabweichung bei " << n << " durchläufen beträgt: " << standardAbweichung << "" << endl;

}

#endif //AWR_PRAKTIKUM_MATHTOOLS_H
