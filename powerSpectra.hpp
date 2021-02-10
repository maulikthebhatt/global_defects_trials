#ifndef POWERSPECTRA_HPP
#define POWERSPECTRA_HPP


void output_powerSpectrum(Field<Imag> & field,
                          string filename,
                          const int nbins,
                          const double physicalSize,
                          const bool deconv = true,
                          const bool dimensionless = false,
                          const bool normFT = false,
                          const bool factor_norm = false)
{

  if(field.nMatrix()!=1)
  {
    COUT<< "output_powerSpectrum not suported for rank 3 tensor.\n Aborting to write: "<<filename<<endl;
    return;
  }

  fstream file;

  if(parallel.rank()==0)
  {
    file.open( filename.c_str() , std::fstream::out | std::fstream::trunc);
    if(!file.is_open())
    {
      cout << "cannot open file: " << filename << ", aborting output_powerSpectrum" << endl;
    }
  }
  
  rKSite k(field.lattice());
  long npts = field.lattice().size(1);

  double * k_arr     = new double[nbins];
  double * kstd_arr  = new double[nbins];
  double * ps_arr    = new double[nbins];
  double * psstd_arr = new double[nbins];
  int * n_arr        = new int[nbins];
  double * s_arr     = new double[npts];
  double * k2_arr    = new double[npts];

  double normFT_factor;
  if(normFT)
  {
    normFT_factor = npts*npts*npts*npts*npts*npts;
  }
  else
  {
    normFT_factor = 1.0;
  }

  double norm_factor;
  if(factor_norm)
  {
    norm_factor = pow(physicalSize,3);
  }
  else
  {
    norm_factor = 1.0;
  }

  double w;
  double k2, k2max,sqrtk2;
  double s2;
  Imag ps;

  double kfactor = 2.0*M_PI/physicalSize;

  int i;



  for (i = 0; i <= npts/2; i++)
  {
    k2_arr[i] = 2. * M_PI * (double)i;
    k2_arr[i] /= physicalSize;
    k2_arr[i] *= k2_arr[i];
  }
  for (; i<npts; i++)
  {
    k2_arr[i] = 2. * M_PI * (double)(npts-i);
    k2_arr[i] /= physicalSize;
    k2_arr[i] *= k2_arr[i];
  }
  k2max = 3.0 * k2_arr[npts/2];

  s_arr[0] = 1.0;
  if(deconv)
	{
		for (i = 1; i <= npts / 2; i++)
		  s_arr[i] = sin(M_PI * (double)i / (double)npts) * (double)npts / (M_PI * (double)i);
	}
	else
	{
		for (i = 1; i <= npts / 2; i++) s_arr[i] = 1.0;
	}
	for (; i < npts; i++) s_arr[i] = s_arr[npts-i];

  for (i = 0; i < nbins; i++)
	{
		k_arr[i]     = 0.0;
		kstd_arr[i]  = 0.0;
		ps_arr[i]    = 0.0;
		psstd_arr[i] = 0.0;
		n_arr[i]     = 0;
	}

  for (k.first(); k.test(); k.next())
	{
		if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0) continue;
		else if(k.coord(0) == 0 || k.coord(0) == npts/2)
			w = 1.0;
		else
			w = 2.0;

		k2 = k2_arr[k.coord(0)] + k2_arr[k.coord(1)] + k2_arr[k.coord(2)];
		s2 = s_arr[k.coord(0)] * s_arr[k.coord(1)] * s_arr[k.coord(2)];
		s2 *= s2;

    if (field.symmetry() == LATfield2::symmetric)
		{
			ps = field(k, 0, 1).norm();
			ps += field(k, 0, 2).norm();
			ps += field(k, 1, 2).norm();
			ps *= 2.0;
			ps += field(k, 0, 0).norm();
			ps += field(k, 1, 1).norm();
			ps += field(k, 2, 2).norm();
		}
		else
		{
			ps = Imag(0., 0.);
			for (i = 0; i < field.components(); i++)
			   ps += field(k, i).norm();
		}


    sqrtk2 = std::sqrt(k2);
    i = std::floor( (double)nbins * std::sqrt(k2 / k2max) );
		if (i<nbins)
		{
			k_arr[i]     += w * sqrtk2;
			kstd_arr[i]  += w * k2;
      if(dimensionless)
      {
          ps_arr[i]    += w * ps.real() / (sqrtk2*k2)  / s2 / normFT_factor * norm_factor;
          psstd_arr[i] += w * ps.real() * ps.real() / (k2*k2*k2) / (s2*s2) / normFT_factor / normFT_factor * norm_factor * norm_factor;
      }
      else
      {
          ps_arr[i]    += w * ps.real() / s2 / normFT_factor * norm_factor;
          psstd_arr[i] += w * ps.real() * ps.real() / (s2*s2) / normFT_factor / normFT_factor * norm_factor * norm_factor;
      }
      n_arr[i]     += w;
		}

  }

  parallel.sum(k_arr, nbins);
	parallel.sum(kstd_arr, nbins);
	parallel.sum(ps_arr, nbins);
	parallel.sum(psstd_arr, nbins);
	parallel.sum(n_arr, nbins);

	for (i = 0; i < nbins; i++)
	{
		if (n_arr[i] > 0)
		{
      k_arr[i] /= (double)n_arr[i];
      kstd_arr[i] = sqrt(kstd_arr[i] / (double)n_arr[i] - k_arr[i] * k_arr[i]);
			if (!isfinite(kstd_arr[i])) kstd_arr[i] = 0.;
      ps_arr[i] /= (double)n_arr[i];
      psstd_arr[i] = sqrt(psstd_arr[i] / (double)n_arr[i] - ps_arr[i] * ps_arr[i]);
			if (!isfinite(psstd_arr[i])) psstd_arr[i] = 0.;
		}
	}

  if(parallel.rank()==0)
  {
    file<<"npts: "<< npts
        << ", physicalSize: "<<physicalSize
        <<endl;


    for (i = 0; i < nbins; i++)
  	{
      file<<i<<"   "<< k_arr[i]<< "   "
          << kstd_arr[i] << "   "
          << ps_arr[i] << "   "
          << psstd_arr[i] << "   "
          << n_arr[i] <<endl;
    }
    file.close();
  }


  delete[] k_arr;
  delete[] kstd_arr;
  delete[] ps_arr;
  delete[] psstd_arr;
  delete[] n_arr;
  delete[] s_arr;
  delete[] k2_arr;


}

#endif
