//---------------------------------------------------------------------------

#include "stdafx.h"

#include "CostFunction.h"

//---------------------------------------------------------------------------

CostFunction::CostFunction(const CImg<unsigned char>& target, const CImg<unsigned char>& source)
{
	this->target = target;
	this->source = source;
}

// p = paramètres de la transformation rigide 
// p[0] = translation en x
// p[1] = translation en y
// p[2] = angle de rotation (en degrés)
double ChamferDistance::Evaluate(double* p)
{
	double angle = std::acos(-1.0) * p[st + 2] / 180;
	double cosa = std::cos(angle), sina = std::sin(angle);
	
	// Pour tous les points du contour de l'image source
	// appliquer la transformation rigide 
	// et interpoler la carte de distance au point transformé
	double sum = 0.0;
	cimg_forXY(source, x, y)
	{
		if (source(x, y))
		{
			// Passage dans le repère du centre de l'image source
			double xcs = static_cast<double>(x) - source.dimx() / 2;
			double ycs = static_cast<double>(y) - source.dimy() / 2;

			// Rotation autour du centre de l'image source d'angle p[2]
			double xct = xcs * cosa + ycs * sina;
			double yct = -xcs * sina + ycs * cosa;

			// Translation
			xct += p[st + 0];
			yct += p[st + 1];

			// Passage dans le repère de l'image cible
			double xt = xct + target.dimx() / 2;
			double yt = yct + target.dimy() / 2;

			// Interpolation bicubique de l'image cible en (xt, yt)
			if (xt >= 0 && xt < target.dimx() && yt >= 0 && yt < target.dimy())
				sum += target.cubic_atXY(xt, yt);
			else
			{
				unsigned char min;
				sum += target.maxmin(min);
			}
		}
	}

	return sum;
}

MutualInformation::MutualInformation(const CImg<unsigned char>& target, const CImg<unsigned char>& source, int bin, int subsample):CostFunction(target, source)
{
	// Précalcul de l'histogramme de l'image à réorienter
	this->bin = bin;
	this->subsample = subsample;
}

double MutualInformation::Evaluate(double* p)
{
	double angle = std::acos(-1.0) * p[st + 2] / 180;
	double cosa = std::cos(angle), sina = std::sin(angle);

	CImg<double> h_joint(bin, bin), h_source(bin), h_target(bin);
	h_joint.fill(0.0);
	h_source.fill(0.0);
	h_target.fill(0.0);

	// Parcours des pixels de l'image source avec un pas subsample
	int np = 0;
	for (int x = 0; x < source.dimx(); x += subsample)
		for (int y = 0; y < source.dimy(); y += subsample)
		{
			// Application de la transformation rigide à chaque pixel
			// Passage dans le repère du centre de l'image source
			double xcs = static_cast<double>(x) - source.dimx() / 2;
			double ycs = static_cast<double>(y) - source.dimy() / 2;

			// Rotation autour du centre de l'image source d'angle p[2]
			double xct = xcs * cosa + ycs * sina;
			double yct = -xcs * sina + ycs * cosa;

			// Translation
			xct += p[st + 0];
			yct += p[st + 1];

			// Passage dans le repère de l'image cible
			double xt = xct + target.dimx() / 2;
			double yt = yct + target.dimy() / 2;
			
			if (xt >= 0 && xt < target.dimx() && yt >= 0 && yt < target.dimy())
			{
				// Niveau de gris de l'image source
				unsigned char ng_source = static_cast<unsigned char>(source(x, y) * (bin - 1) / 255);

				// Niveau de gris de l'image cible
				unsigned char ng_target = static_cast<unsigned char>(target.cubic_atXY(xt, yt) * (bin - 1) / 255);
				
				// Mise à jour de l'histogramme joint
				h_joint(ng_source, ng_target)++;

				// Mise à jour de l'histogramme marginal de l'image source
				h_source(ng_source)++;

				// Mise à jour de l'histogramme marginal de l'image cible
				h_target(ng_target)++;

				np++;
			}
		}

	h_joint /= np;
	h_source /= np;
	h_target /= np;

	// Calcul de l'information mutuelle
	double MI = 0.0;
	for (int a = 0; a < bin; a++)
		for (int b = 0; b < bin; b++)
			if (h_joint(a, b) != 0.0 && h_source(a) != 0.0 && h_target(b) != 0.0)
				MI += h_joint(a, b) * std::log(h_joint(a, b) / h_source(a) / h_target(b));

	return -MI;
}