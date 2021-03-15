
#include "../headers/GlFibreCorrection.h"


/**
 * @brief Normalizes fiber and sheet directions of a <CODE>.lon</CODE> file so that they have norm 1.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.lon</CODE> file to normalize.
 * @param fileName (<CODE>string</CODE>) name of the <CODE>.lon</CODE> file to normalize.
 */
void normalizeLon(const std::string &path,const std::string &fileName)
{
	std::string outPath;
	/* Read .lon file */
	lon fibres = ReadLon( path, fileName );

	lon fibresNorm = fibres;

	int kF = 0;
	int kS = 0;

	/* Normalize fibres */
	for( int i = 0; i < fibres.x1.size(); i++ )
	{
		double normF = sqrt( pow( fibres.x1[i], 2 ) + pow( fibres.y1[i], 2 ) + pow( fibres.z1[i], 2 ) );
		if( (normF >= 1.0 + 1e-5) || (normF <= 1.0 - 1e-5) )
		{
			fibresNorm.x1[i] = fibres.x1[i]/normF;
			fibresNorm.y1[i] = fibres.y1[i]/normF;
			fibresNorm.z1[i] = fibres.z1[i]/normF;
			kF = kF + 1;
		}
		double normS = sqrt( pow( fibres.x2[i], 2 ) + pow( fibres.y2[i], 2 ) + pow( fibres.z2[i], 2 ) );
		if( (normS >= 1.0 + 1e-5) || (normS <= 1.0 - 1e-5) )
		{
			fibresNorm.x2[i] = fibres.x2[i]/normS;
			fibresNorm.y2[i] = fibres.y2[i]/normS;
			fibresNorm.z2[i] = fibres.z2[i]/normS;
			kS = kS + 1;
		}
	}

	std::fstream fOutLon;
	outPath = path + "/" + fileName + "_norm" + ".lon";
	write2file(outPath, fOutLon );
	fOutLon << 2 << std::endl;

	for( int i = 0; i < fibres.x1.size(); i++ )
	{
		fOutLon << fibresNorm.x1[i] << " " << fibresNorm.y1[i] << " " << fibresNorm.z1[i] << 
		    " " << fibresNorm.x2[i] << " " << fibresNorm.y2[i] << " " << fibresNorm.z2[i] << std::endl;
	}

	fOutLon.close();

	std::cout << "Corrected " << kF << " fibre vectors" << std::endl;
	std::cout << "Corrected " << kS << " sheet vectors" << std::endl;

}


std::vector<double> matVec_prod(const std::vector<std::vector<double>> &M,const std::vector<double> &v )
{	
	std::vector<double> p(3);
	double product = 0.0;

	for (int i = 0; i < 3; i++)
	{
	    for (int j = 0; j < 3; j++)
	    {
	        product += M[i][j] * v[j];
	        p[i] = product;
	    }
	    product = 0.0;
	}
	return p;

}

double dot_prod(const std::vector<double> &u,const std::vector<double> &v )
{	
	double p = 0.0;

	for (int i = 0; i < 3; i++)
	{
		p = p + u[i]*v[i];
	}
	return p;
}

std::vector<double> cross_prod(const std::vector<double> &u,const std::vector<double> &v )
{	
	std::vector<double> p(3);

    p[0] = u[1] * v[2] - u[2] * v[1]; 
    p[1] = u[2] * v[0] - u[0] * v[2]; 
    p[2] = u[0] * v[1] - u[1] * v[0]; 

    return p;
}

std::vector<double> scalar_prod(const std::vector<double> &u,const double &k )
{	
	std::vector<double> p(3);

	for (int i = 0; i < 3; i++)
	{
		p[i] = u[i] * k;
	}
	return p;
}

std::vector<double> orthogonalise(const std::vector<double> &f,const std::vector<double> &s )
{	

	/* cross product */
	std::vector<double> c(3), axis(3);
	double norm_c, theta, d;
	double pi = 3.14159265358979323846;

	c = cross_prod( f, s );
	d = dot_prod( f, s );

	norm_c = sqrt( dot_prod( c, c ) ); 
	axis = scalar_prod( c, 1/norm_c );
	
	theta = pi/2 - std::acos( d );

	/* Define the rotation matrix */
	std::vector<std::vector<double>> R(4, std::vector<double>(4));
	
	R[0][0] = pow(axis[0],2) + std::cos(theta) * (1 - pow(axis[0],2));
	R[0][1] = (1 - std::cos(theta)) * axis[0] * axis[1] - axis[2] * std::sin(theta);
	R[0][2] = (1 - std::cos(theta)) * axis[0] * axis[2] + axis[1] * std::sin(theta);
	R[1][0] = (1 - std::cos(theta)) * axis[0] * axis[1] + axis[2] * std::sin(theta);
	R[1][1] = pow(axis[1],2) + std::cos(theta) * (1 - pow(axis[1],2));
	R[1][2] = ( 1 - cos(theta)) * axis[1] * axis[2] - axis[0] * std::sin(theta);
	R[2][0] = ( 1 - cos(theta)) * axis[0] * axis[2] - axis[1] * std::sin(theta);
	R[2][1] = ( 1 - cos(theta)) * axis[1] * axis[2] + axis[0] * std::sin(theta);
	R[2][2] = pow(axis[2],2) + std::cos(theta) * (1 - pow(axis[2],2));

	/* Rotate sheet vector */
	std::vector<double> s_corrected(3);	
	s_corrected = matVec_prod( R, s );

	return s_corrected;
}

void orthogonalise_lon(const std::string &mshPath,const std::string &mshName )
{

	/* Read file */
	lon fibres = ReadLon( mshPath, mshName );

	lon fibres_corrected = fibres;
	std::vector<double> f(3), s(3), s_corrected(3);
	double d;
	int count = 0;

	// #pragma omp parallel for
	for(int i = 0; i < fibres.x1.size(); i++ )
	{
		f[0] = fibres.x1[i]; f[1] = fibres.y1[i]; f[2] = fibres.z1[i];
		s[0] = fibres.x2[i]; s[1] = fibres.y2[i]; s[2] = fibres.z2[i];

		d = dot_prod( f, s );

		if( std::abs(d) >= 1e-05 )
		{
			s_corrected = orthogonalise( f, s );
			fibres_corrected.x2[i] = s_corrected[0];
			fibres_corrected.y2[i] = s_corrected[1];
			fibres_corrected.z2[i] = s_corrected[2];
			count = count+1;
		}
	}

	std::cout << "Corrected sheet direction for " << count << " elements" << std::endl;

	/* Write output */
	std::fstream fOutLon;
    std::string outName = mshPath + "/" + mshName + "_orto" + ".lon";
	write2file( outName, fOutLon );
	fOutLon << 2 << std::endl;

	for( int i = 0; i < fibres.x1.size(); i++ )
	{
		fOutLon << fibres_corrected.x1[i] << " " << fibres_corrected.y1[i] << " " << fibres_corrected.z1[i] << 
		    " " << fibres_corrected.x2[i] << " " << fibres_corrected.y2[i] << " " << fibres_corrected.z2[i] << std::endl;
	}

	fOutLon.close();

}



void GlFibreCorrection(const std::string &mshPath,const std::string &mshName)
{
	
	/* Read mesh in */
	points pts = ReadPts(mshPath,mshName);
	elem tets = ReadElem(mshPath,mshName);
	lon fibres = ReadLon(mshPath,mshName);
	fibres.nDir = 2;

	/* Identify tets with wrong fibres */
	std::vector<int> tets_ind;

	for(int i = 0; i < fibres.x1.size(); i++)
	{
		if(std::abs(fibres.z1[i]) < 1e-06)
		{
			tets_ind.push_back(i);
		}
	}
	std::sort(tets_ind.begin(),tets_ind.end());

	std::cout << "Found " << tets_ind.size() << " to correct" <<std::endl;

	/* Correct fibre direction */
	lon fibres_corrected = fibres;

	#pragma omp parallel for
	for(int i = 0; i < tets_ind.size() ; i++)
	{	
		double r = 2500;   // micrometre
		std::vector<double> g(3), f(3), s(3);
		std::vector<int> sphere, sphere_tets;
		/* Compute centre of gravity of the element */
		g[0] = (pts.x[tets.v1[tets_ind[i]]]+pts.x[tets.v2[tets_ind[i]]]+pts.x[tets.v3[tets_ind[i]]]+pts.x[tets.v4[tets_ind[i]]])/4;
		g[1] = (pts.y[tets.v1[tets_ind[i]]]+pts.y[tets.v2[tets_ind[i]]]+pts.y[tets.v3[tets_ind[i]]]+pts.y[tets.v4[tets_ind[i]]])/4;
		g[2] = (pts.z[tets.v1[tets_ind[i]]]+pts.z[tets.v2[tets_ind[i]]]+pts.z[tets.v3[tets_ind[i]]]+pts.z[tets.v4[tets_ind[i]]])/4;

		/* Find neighbours nodes with a sphere  centred in the centre of gravity of the element */
		for( int j = 0; j < pts.x.size(); j++ )
		{	
			if( (pow(pts.x[j]-g[0],2)+pow(pts.y[j]-g[1],2)+pow(pts.z[j]-g[2],2)) <= pow(r,2) ) 
			{
				sphere.push_back(j);
			}
		}		
	    std::sort(sphere.begin(),sphere.end());

		/* Find neighbours tets with all the nodes in the sphere and take the average of the f and s directions */
		int count = 0;
		double sum_fx = 0;
		double sum_fy = 0;
		double sum_fz = 0;
		double sum_sx = 0;
		double sum_sy = 0;
		double sum_sz = 0;

		for( int k = 0; k < tets.v1.size(); k++ )
		{	
			if( (std::find(sphere.begin(), sphere.end(), tets.v1[k]) != sphere.end()) && 
				(std::find(sphere.begin(), sphere.end(), tets.v2[k]) != sphere.end()) && 
				(std::find(sphere.begin(), sphere.end(), tets.v3[k]) != sphere.end()) && 
				(std::find(sphere.begin(), sphere.end(), tets.v4[k]) != sphere.end()) && 
				(std::find(tets_ind.begin(), tets_ind.end(), k) == tets_ind.end()))
			{
				sphere_tets.push_back(k);
				count = count+1;
				sum_fx = sum_fx + fibres.x1[k];
				sum_fy = sum_fy + fibres.y1[k];
				sum_fz = sum_fz + fibres.z1[k];
				sum_sx = sum_sx + fibres.x2[k];
				sum_sy = sum_sy + fibres.y2[k];
				sum_sz = sum_sz + fibres.z2[k];
			}
		}

		if(count == 0)
		{
			std::cout << "Increase sphere radius..." << std::endl;
		}

		f[0] = sum_fx/count;
		f[1] = sum_fy/count;
		f[2] = sum_fz/count;

		double normF = sqrt(pow(f[0],2)+pow(f[1],2)+pow(f[2],2));
		fibres_corrected.x1[tets_ind[i]] = f[0]/normF;
		fibres_corrected.y1[tets_ind[i]] = f[1]/normF;
		fibres_corrected.z1[tets_ind[i]] = f[2]/normF;

		s[0] = sum_sx/count;
		s[1] = sum_sy/count;
		s[2] = sum_sz/count;


		/* Orthogonalise */
		double d = dot_prod( f, s );
		std::vector<double> s_corrected(3);

		if( std::abs(d) >= 1e-05 )
			s_corrected = orthogonalise( f, s );
		else
        	s_corrected = s;

		double normS = sqrt(pow(s[0],2)+pow(s[1],2)+pow(s[2],2));
		fibres_corrected.x2[tets_ind[i]] = s[0]/normS;
		fibres_corrected.y2[tets_ind[i]] = s[1]/normS;
		fibres_corrected.z2[tets_ind[i]] = s[2]/normS;

		// std::cout << fibres.x1[tets_ind[i]] << " " << fibres.y1[tets_ind[i]] << " " << fibres.z1[tets_ind[i]] << std::endl;
		// std::cout << fibres_corrected.x1[tets_ind[i]] << " " << fibres_corrected.y1[tets_ind[i]] << " " << fibres_corrected.z1[tets_ind[i]] << std::endl;

	}

    std::string outName = mshPath + "/" + mshName + "_corrected" + ".lon";
	std::fstream fOutLon;
	write2file(outName , fOutLon );
	fOutLon << 2 << std::endl;

	for( int i = 0; i < fibres.x1.size(); i++ )
	{
		fOutLon << fibres_corrected.x1[i] << " " << fibres_corrected.y1[i] << " " << fibres_corrected.z1[i] << 
		    " " << fibres_corrected.x2[i] << " " << fibres_corrected.y2[i] << " " << fibres_corrected.z2[i] << std::endl;
	}

	fOutLon.close();
}

