using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coreflood
{

	public class RCore
	{
		public readonly double a;
		public readonly double b;
		public double c { get; set; }
		public readonly double porosity;
		public readonly double diameter;
		public readonly double length;
		public readonly double permeability;
		public readonly double krwMAX;
		public readonly double kroMAX;
		public readonly double pcMAX;
		public readonly double swi;
		public readonly double sor;
		public readonly double Pct;
		public readonly static double pi = 3.14159265358979323846; //const in c++
		public readonly string wettability;

		public RCore(CoreArgs cr, string _w = "W") {
			porosity = cr.porosity;
			diameter = cr.diameter; 
			length = cr.length; 
			permeability = cr.permeability; 
			wettability = _w; swi = cr.swi; 
			sor = cr.sor; 
			a = cr.exp_krw; 
			b = cr.exp_kro;
			c = cr.exp_pc; 
			krwMAX = cr.krwMax; 
			kroMAX = cr.kroMax; 
			pcMAX = cr.pcMax;
			Pct = cr.pct;
		}
		public double area() => 0.25 * pi * diameter * diameter;
		public double volume() => area() * length;
		public double volPoroso() => volume() * porosity;
		public double krw(double s) {
			double r;
			if (s >= 1 - sor) { r = krwMAX; }
			else if (s <= swi) { r = 0; }
			else {
				r = krwMAX * Math.Pow(((s - swi) / (1 - sor - swi)), a);
			}
			return r;	
		}
		public double kro(double s) {
			double r;
			if (s >= 1 - sor) { r = 0; }
			else if (s <= swi) { r = kroMAX; }
			else
			{
				r = kroMAX * Math.Pow(((1 - sor - s) / (1 - sor - swi)), a);
			}
			return r;
		}
		public double lambda_w(double s, double uw) => krw(s) / uw;
		public double lambda_o(double s, double uo) => kro(s) / uo;
		public double lambda_t(double s, double uw, double uo) => lambda_w(s, uw) + lambda_o(s, uo);
		public double pc(double s){
			double r;
			if (s >= 1 - sor) { r = Pct; }
			else
			{
				r = pcMAX * Math.Pow(((1 - sor - s) / (1 - sor - swi)), c) + Pct;
			}
			return r;
		}
		public double fw(double s, double uw, double uo) => lambda_w(s, uw) / lambda_t(s, uw, uo);

		public void PlotKr() { }
		public void PlotPc() { }

		//friend class IMPES;
		//friend class IMPES_Parallel;
		//friend class IMPES_ParallelCV;
	}
}
