//  définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "CostFunction.h"
#include "Optimisation.h"
#include "../CImg.h"

using namespace std;
using namespace cimg_library;

CImg<unsigned char> getSkullContours(const CImg<unsigned char>& img, const unsigned char threshold, const int radius)
{
	// Seuillage au niveau threshold
	CImg<unsigned char> bin(img.dimx(), img.dimy());
	bin = img.get_threshold(threshold);

	// Fermeture avec un élément structurant circulaire de rayon radius
	CImg<unsigned char> es(2 * radius + 1, 2 * radius + 1);
	for (int x = -radius; x <= radius; x++)
		for (int y = -radius; y < radius; y++)
		{
			if (x * x + y * y < radius * radius)
				es(x + radius, y + radius) = 1;
			else
				es(x + radius, y + radius) = 0;
		}
	bin.dilate(es);
	bin.erode(es);
	
	// Détection du contour de l'image seuillée
	CImg<unsigned char> skull(img.dimx(),img.dimy());
	cimg_forXY(bin, x, y)
	{
		skull(x, y) =   ( bin(x, y) == 1 ) &&
						( (bin(x + 1, y) == 0) || 
						  (bin(x - 1, y) == 0) || 
						  (bin(x, y + 1) == 0) || 
						  (bin(x, y - 1) == 0) );
	}
	return skull;
}

void ChamferMatching(const CImg<unsigned char>& A, const double threshold_rateA, const CImg<unsigned char>& B, double threshold_rateB, CImg<double>& translation, double& rotation, int optim)
{
	// Calcul de la distance du chanfrein à l'image A
	unsigned char maxA = 0;
	cimg_forXY(A,x,y) {if (maxA<A(x,y)) maxA = A(x,y);}
	unsigned char thresholdA = (unsigned char)(threshold_rateA * maxA);
	CImg<unsigned char> skullA = getSkullContours(A,thresholdA,9);
	skullA.distance(1);
	//skullA.display();

	// Calcul des contours du crâne dans l'image B
	unsigned char maxB = 0;
	cimg_forXY(B,x,y) {if (maxB<B(x,y)) maxB = B(x,y);}
	unsigned char thresholdB = (unsigned char)(threshold_rateB * maxB);
	CImg<unsigned char> skullB = getSkullContours(B,thresholdB,9);
	skullB.display();

	ChamferDistance CD(skullA,skullB);

	CImg<double> res;
	if (optim==0) res = QuasiNewton( CD, translation, rotation );
	if (optim==1) res = StochasticClustering( CD, translation, rotation, 20, 10 );
	if (optim==2) res = SimulatedAnnealing( CD, translation, rotation );
}

void MutualInformationMatching(const CImg<unsigned char>& A, const CImg<unsigned char>& B, CImg<double>& translation, double& rotation, int bin, int subsample, bool optim)
{
	MutualInformation MI(A,B,bin,subsample);

	CImg<double> res;
	if (optim==0) res = QuasiNewton( MI, translation, rotation );
	if (optim==1) res = StochasticClustering( MI, translation, rotation, 20, 10 );
	if (optim==2) res = SimulatedAnnealing( MI, translation, rotation );
}

int _tmain(int argc, _TCHAR* argv[])
{
 // Ouverture des images d'entrée
 CImg<unsigned char> A = CImg<>("../Img/MRI.bmp").channel(0);
 CImg<unsigned char> B0 = CImg<>("../Img/MRI.bmp").channel(0);
 //CImg<unsigned char> A = CImg<>("../Img/MRI.bmp").channel(0);
 //CImg<unsigned char> B0 = CImg<>("../Img/TEMP.bmp").channel(0);
 CImg<unsigned char> B(B0);
 
 CImg<double> translation(2);
 double rotation;

 cout << "Translation en x : ";
 cin >> translation(0);
 cout << "Translation en y : ";
 cin >> translation(1);
 cout << "Rotation en degres : ";
 cin >> rotation;
 
 // Simulation d'une transformation rigide sur l'image cible
 // si oui, les valeurs initiales de paramètres sont mis à zéro
 // si non, les valeurs initiales sont les valeurs choisies
 int choice;
 cout << "Simulation d'une transformation rigide (0:non, 1:oui) : ";
 cin >> choice;
 if (choice)
 {
	A.rotate(-rotation,2,2);
	A.translate(translation(0),translation(1)); 
	translation.fill(0.);
	rotation = 0.;
 }

 char str[256];

 cout << "Choix du critere de recalage (0:contour, 1:intensite) : ";
 cin >> choice;
 int optim;
 cout << "Choix de la methode d'optimisation (0:Quasi-Newton, 1:Stochastic Clustering, 2:simulated annealing) : ";
 cin >> optim;

 // Recalage optimisant la distance du chanfrein du contour
 if (choice == 0)
 {
	strcpy(str,"Distance du chanfrein");
	ChamferMatching(A,.35,B0,.35,translation,rotation,optim);
 }
 
 // Recalage utilisant l'information mutuelle comme mesure de similarité
 if (choice == 1)
 {
	strcpy(str,"Information mutuelle");
	int bin;
	cout << "Nombre de bins pour le calcul des histogrammes : ";
	cin >> bin;
	int subsample;
	cout << "Facteur de sous-echantillonnage des images : ";
	cin >> subsample;
	MutualInformationMatching(A,B0,translation,rotation,bin,subsample,optim);
 }

 cout << "Translation estimee en x : " << translation(0) << endl;
 cout << "Translation estimee en y : " << translation(1) << endl;;
 cout << "Rotation en degres : " << rotation << endl;

 B.rotate(-rotation,2,2);
 B.translate(translation(0),translation(1)); 

 A.resize(512,512);
 CImgDisplay disp_A(A,"Image cible");
 B0.resize(512,512);
 CImgDisplay disp_B0(B0,"Position initiale");
 B.resize(512,512);
 CImgDisplay disp_B(B,"Position finale",1);

 // Affichages superposés en damier
 int min_dimx = (A.dimx()<B.dimx())?A.dimx():B.dimx();
 int min_dimy = (A.dimy()<B.dimy())?A.dimy():B.dimy();
 int offsetA_x = .5*(A.dimx()-min_dimx), offsetA_y = .5*(A.dimy()-min_dimy);
 int offsetB_x = .5*(B.dimx()-min_dimx), offsetB_y = .5*(B.dimy()-min_dimy);
 A.crop(offsetA_x,offsetA_y,A.dimx()-offsetA_x-1,A.dimy()-offsetA_y-1);
 B.crop(offsetB_x,offsetB_y,B.dimx()-offsetB_x-1,B.dimy()-offsetB_y-1);


 CImg<unsigned char> AB(min_dimx, min_dimy);
 int nx = 8, ny = 8, dx = min_dimx / nx, dy = min_dimy / ny;

 for (int i = 0; i < min_dimx; i++)
 {
	for (int j = 0; j < min_dimy; j++)
	{
		// Case image A
		if ( ((i/dx)%2) ^ ((j/dy)%2) )
		{
			AB(i,j) = A(i,j);
		}
		// Case image B
		else
		{
			AB(i,j) = B(i,j);
		}
	}
 }

 CImgDisplay disp_AB(AB,"Affichage superpose",1);

 while (!disp_A.is_closed ) disp_A.wait();

 return 0;
}
