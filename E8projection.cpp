#include "E8projection.h"

/*
	Projected distances are divided by 2 to avoid problems with .cif created in Materials Studio (see line 198, 215)
*/
E8projection::E8projection(string filename)
{
	/*input parser skips blank lines
	and lines beginning with #
	*/
	E8D3D(3,8);
	current3Dpoint(3);
	current8Dpoint(8);
	shiftvector(8);
	inputfile.open(filename);
	string temp;
	if(inputfile.is_open())
	{
		while(!inputfile.eof())
		{
			getline(inputfile,temp);
			if(temp[0]=='#')
			{
				/* do nothing */;
			}
			if(temp.compare("LAYERS")==0)
			{
				getline(inputfile,temp);
				sscanf(temp.data(), "%i", &LAYERS);
			}
			else if(temp.compare("CUTOFF")==0)
			{
				getline(inputfile,temp);
				sscanf(temp.data(), "%f", &MAXdist);
			}
			else if(temp.compare("PROJECTIONMATRIX")==0)
			{
				float f[3][8];
				for(int i=0; i < 3; i++)
				{
					getline(inputfile,temp);
					sscanf(temp.data(), "%f %f %f %f %f %f %f %f", &f[i][0],&f[i][1],&f[i][2],&f[i][3],&f[i][4],&f[i][5],&f[i][6],&f[i][7]);						
				}
					E8D3D << f[0][0] << f[0][1] << f[0][2] << f[0][3] << f[0][4] << f[0][5] << f[0][6] << f[0][7] << endr
						<< f[1][0] << f[1][1] << f[1][2] << f[1][3] << f[1][4] << f[1][5] << f[1][6] << f[1][7] << endr
						<< f[2][0] << f[2][1] << f[2][2] << f[2][3] << f[2][4] << f[2][5] << f[2][6] << f[2][7] << endr;
			}
			else if(temp.compare("SHIFTVECTOR")==0)
			{
				float f[8];
				getline(inputfile,temp);
				sscanf(temp.data(), "%f %f %f %f %f %f %f %f", &f[0],&f[1],&f[2],&f[3],&f[4],&f[5],&f[6],&f[7]);
				shiftvector << f[0] << f[1] << f[2] << f[3] << f[4] << f[5] << f[6] << f[7]; 
			}
			else
				/* do nothing */;
		}
	}
	else
	{
		cout << "Error opening " << filename << endl;
		exit(-1);
	}
	E8D3Dtr = trans(E8D3D);
	basis3D = E8D3D * E8D3Dtr;
	scale3D = inv(basis3D);
	inputfile.close();
}
E8projection::E8projection(float t)
{
	E8D3D(3,8);
	E8D3D << 0.00 << 0.00 << -2*t << 0.00 << 1.00 << 1.00 << -1.00 << -1.00 << endr
		<< 0.00 << 0.00 << 0.00 << -2*t << 1.00 << -1.00 << 1.00 << -1.00 << endr
		<< 0.00 << -2*t << 0.00 << 0.00 << 1.00 << -1.00 << -1.00 << 1.00 << endr;

	E8D3Dtr = trans(E8D3D);
	basis3D = E8D3D * E8D3Dtr;
	scale3D = inv(basis3D);
	MAXdist = 0.5;
	LAYERS = 4;
	current3Dpoint(3);
	current8Dpoint(8);
	shiftvector(8);
	shiftvector << 0.36 << 0.00 << 0.00 << 0.00 << 0.11 << 0.11 << 0.11 << 0.11;
}

mat E8projection::orthonormalize(const mat &toorth)
{
	mat orthogonalized(toorth);
	colvec u(3); colvec v(3);
	for(int j = 0; j < (int)orthogonalized.n_cols; j++)
	{
		v = orthogonalized.col(j);
		for(int i = 0; i < (j - 1); i++)
		{
			u = orthogonalized.col(i);
			v = v - (as_scalar(trans(v)*u)/as_scalar(trans(u)*u))*v;
		}
		v = v/norm(v,2);
		for(int i = 0; i < (int)v.n_rows; i++)
		{
			orthogonalized(i,j)=v(i);
		}
	}
	return orthogonalized;
}

void E8projection::project(colvec point8D)
{
	current8Dpoint = point8D-shiftvector;
	current3Dpoint = E8D3D * current8Dpoint;
	current3Dpointcc = scale3D * current3Dpoint;
	currentprojdist = norm(current8Dpoint - E8D3Dtr*current3Dpointcc,2);
}

void E8projection::printbasis3D(string name)
{
	basis3D.print(name);
}

void E8projection::printbasis3D()
{
	basis3D.print();
}

void E8projection::printprojectionmatrix()
{
	E8D3D.print();
}
void E8projection::printprojectionmatrix(string name)
{
	E8D3D.print(name);
}

void E8projection::printshiftvector(string name)
{
	shiftvector.print_trans(name);
}

E8projection::~E8projection(void)
{

}

int E8projection::execute(string output)
{
	ofstream out;
	ofstream outres;
	out.open(output);
	outres.open(output.append(".res"));
	out.precision(5);
	outres.precision(5);
	outres.setf(ios::fixed);
	//Prepare output file
	out << __DATE__ << " " << __TIME__ << endl;
	out << "PROJECTION MATRIX" << endl << E8D3D << endl;
	out << "ORIGIN SHIFT" << endl << trans(shiftvector) << endl;
	out << "LAYERS" << endl << LAYERS << endl;
	out << "CUTOFF" << endl << MAXdist << endl;
	out << "COODINATES" << endl;
	//Prepare SHELX file
	outres << "TITL " << output << endl;
	outres << "CELL " << "1.5478 " << norm(basis3D.row(0),2) << " " << norm(basis3D.row(1),2) << " " << norm(basis3D.row(2),2) << " "
		<< std::acos((basis3D.row(1) * trans(basis3D.row(2))/(norm(basis3D.row(1),2)*norm(basis3D.row(2),2))).at(0,0))*(180/math::pi()) << " " 
		<< std::acos((basis3D.row(0) * trans(basis3D.row(2))/(norm(basis3D.row(0),2)*norm(basis3D.row(2),2))).at(0,0))*(180/math::pi()) << " " 
		<< std::acos((basis3D.row(0) * trans(basis3D.row(1))/(norm(basis3D.row(0),2)*norm(basis3D.row(1),2))).at(0,0))*(180/math::pi()) << endl; 
	outres << "LATT -1" << endl;
	outres << "SFAC Cu" << endl;
	colvec vec(8);
	colvec halfshift(8);
	halfshift << 0.5 << 0.5 << 0.5 << 0.5 << 0.5 << 0.5 << 0.5 << 0.5;
	int totalpoints=0;
	int i[8];
	for(i[0]=-LAYERS; i[0] <= LAYERS; i[0]++){
		for(i[1]=-LAYERS; i[1] <= LAYERS; i[1]++){
			for(i[2]=-LAYERS; i[2] <= LAYERS; i[2]++){
				for(i[3]=-LAYERS; i[3] <= LAYERS; i[3]++){
					for(i[4]=-LAYERS; i[4] <= LAYERS; i[4]++){
						for(i[5]=-LAYERS; i[5] <= LAYERS; i[5]++){
							for(i[6]=-LAYERS; i[6] <= LAYERS; i[6]++){
								for(i[7]=-LAYERS; i[7] <= LAYERS; i[7]++){
									if((i[0] + i[1] + i[2] + i[3] + i[4] + i[5] + i[6] + i[7])%2 == 0)
									{
										vec << (float)i[0] << (float)i[1] << (float)i[2] << (float)i[3] << (float)i[4] << (float)i[5] << (float)i[6] << (float)i[7];
										project(vec);
										if(currentprojdist < MAXdist && current3Dpointcc[0]< 1.00 && current3Dpointcc[0]>=0.00 &&
																			current3Dpointcc[1]< 1.00 && current3Dpointcc[1]>=0.00 &&
																			current3Dpointcc[2]< 1.00 && current3Dpointcc[2]>=0.00)
										{	/*
											vec.print("E8 point: ");
											current3Dpoint.print("3D point: ");
											current3Dpointcc.print("3D point CC: ");
											cout << "Projected Distance: " << currentprojdist << endl;
											cout << "*****" << endl;
											*/
											totalpoints++;
											out << totalpoints << endl << trans(vec) << trans(current3Dpoint) << trans(current3Dpointcc) << "\t" << currentprojdist << endl;
											outres << totalpoints<< setw(5) << "1 " << current3Dpointcc.at(0,0) << " " <<  current3Dpointcc.at(1,0) << " " <<  current3Dpointcc.at(2,0) << " " << "11.00000 "<< currentprojdist  << endl;
										}
										vec += halfshift;
										project(vec);
										if(currentprojdist < MAXdist && current3Dpointcc[0]< 1.00 && current3Dpointcc[0]>=0.00 &&
																			current3Dpointcc[1]< 1.00 && current3Dpointcc[1]>=0.00 &&
																			current3Dpointcc[2]< 1.00 && current3Dpointcc[2]>=0.00)
										{
											/*
											vec.print("E8 point: ");
											current3Dpoint.print("3D point: ");
											current3Dpointcc.print("3D point CC: ");
											cout << "Projected Distance: " << currentprojdist << endl;
											cout << "*****" << endl;
											*/
											totalpoints++;
											out << totalpoints << endl << trans(vec) << trans(current3Dpoint) << trans(current3Dpointcc) << "\t" << currentprojdist << endl;
											outres << totalpoints<< setw(5) << "1 " << current3Dpointcc.at(0,0) << " " <<  current3Dpointcc.at(1,0) << " " <<  current3Dpointcc.at(2,0) << " " << "11.00000 "<< currentprojdist  << endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	outres << "END" << endl;
	out.close();
	outres.close();
	return totalpoints;
}
