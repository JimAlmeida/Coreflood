using System;
using System.Collections.Generic;
using System.Linq;

namespace Coreflood
{
	public class IMPES {
		protected bool drenagem;
		protected bool embebicao;
		protected int nbl;
		protected double uw;
		protected double uo;
		protected double Pl;
		protected double Qinj;
		protected double DSMAX;
		protected double h;
		protected double vpi1;
		protected double vpi2;
		protected double t_drenagem;
		protected double t_embebicao;
		protected double tempo_sim;
		protected double si;
		protected RCore core;

		protected List<List<double>> pc;
		protected List<double[]> p;
		protected List<List<double>> s;
		protected List<double> g;

		protected List<string> tempo_plot;
		protected List<string> tempo_plot_emb;
		protected List<List<double>> Plot_s;
		protected List<List<double>> Plot_s_emb;
		protected List<double[]> Plot_p;
		protected List<double[]> Plot_p_emb;
		protected List<double> Plot_h;

		public IMPES(RCore ncore, IMPESArgs ia){
			core = ncore;
			uw = ia.water_viscosity;
			uo = ia.oil_viscosity;
			Pl = ia.atmospheric_pressure;
			Qinj = ia.injection_flow;
			DSMAX = ia.maximum_saturation_step;
			nbl = ia.number_of_blocks;
			vpi1 = ia.pore_volume_injected_in_drainage;
			vpi2 = ia.pore_volume_injected_in_imbibition;
			si = ia.initial_saturation;

			Tuple<double, double> tempos = tempoInj();
			t_drenagem = tempos.Item1;
			t_embebicao = tempos.Item2;

			pc = new();
			p = new();
			s = new();
			g = new();
			tempo_plot = new();
			tempo_plot_emb = new();
			Plot_s = new();
			Plot_p = new();
			Plot_s_emb = new();
			Plot_p_emb = new();
			Plot_h = new();
	}

		protected Tuple<double, double> tempoInj() {
			double param1 = core.volPoroso() / Qinj;
			return new Tuple<double, double>(vpi1 * param1, vpi2 * param1);
		}
	
		public double u() {
			return Qinj / (core.area());
		}
		public double dt() {
			List<double> absG = g.ConvertAll(x => Math.Abs(x));
			double t = DSMAX / absG.Max();

			if (Double.IsInfinity(t)) return DSMAX;
			
			if (drenagem && tempo_sim + t > t_drenagem) return t_drenagem - tempo_sim;
			else if (embebicao && tempo_sim + t > t_embebicao) return t_embebicao - tempo_sim;
			else if (t > 20) return t / 100;
			else return t;
		}

		public double[] thomasSolver(ref double[] a, ref double[] b, ref double[] c, ref double[] d) {
			int N = d.Count();

			c[0] = c[0] / b[0];
			d[0] = d[0] / b[0];

			double m;
			for (int i = 1; i < N; i++)
			{
				m = 1.0 / (b[i] - a[i] * c[i - 1]);
				c[i] = c[i] * m;
				d[i] = (d[i] - a[i] * d[i - 1]) * m;
			}
			for (int i = N - 1; i-- > 0;)
			{
				d[i] = d[i] - (c[i] * d[i + 1]);
			}
			return d;
		}

		public virtual void solverManifold() {
			double[] s0 = new double[nbl];
			Array.Fill(s0, si);
			double[] pc0 = new double[nbl];
			Array.Fill(pc0, core.Pct);

			s.Add(s0.ToList());
			Plot_s.Add(s0.ToList());
			tempo_plot.Add("0");
			pc.Add(pc0.ToList());
			
			int n = 1;
			double t = 0; 
			double t_fim = 0;
			drenagem = true;

			solverPressure(0); //+
			for (double t_sim = 0; t_sim < t_drenagem; t_sim += t)
			{
				t = solverSaturation1(n);
				for (int j = 0; j < nbl; j++)
				{
					pc0[j] = core.pc(s[n][j]);
				}
				pc.Add(pc0.ToList()); //+
				solverPressure(n); //+
				n += 1;

				bool cond1 = n < 3000 && n % 500 == 0;
				bool cond2 = n > 3000 && n % 2000 == 0;
				if (cond1 || cond2)
				{
					Plot_p.Add(p[n - 1]);
					Plot_s.Add(s[n - 1]);
					tempo_plot.Add(t_fim.ToString());
					simLogger(n, t, t_sim);
				}
				t_fim = t_sim;
			}

			drenagem = false;
			//std::cout << "TEMPO FINAL DRENAGEM: " << (clock() - time_drenagem) / CLOCKS_PER_SEC << '\n';

			double[] semb = new double[nbl];
			Array.Fill(semb, core.swi);
			
			s.Add(semb.ToList());
			for (int j = 0; j < nbl; j++)
			{
				pc0[j] = core.pcMAX;
			}

			pc.Add(pc0.ToList());
			solverPressure(n);
			n += 1;
			Plot_s.Add(s[n - 1]);
			Plot_p.Add(p[n - 1]);
			Plot_s_emb.Add(s[n - 1]);
			Plot_p_emb.Add(p[n - 1]);
			tempo_plot_emb.Add(t_drenagem.ToString());
			tempo_plot.Add(t_drenagem.ToString());
			simLogger(n, t, t_fim);

			embebicao = true;
			if (embebicao)
			{
				for (double t_sim = t_drenagem; t_sim < t_embebicao; t_sim += t)
				{
					t = solverSaturation1(n);
					for (int j = 0; j < nbl; j++)
					{
						pc0[j] = core.pc(s[n][j]);
					}
					pc.Add(pc0.ToList());
					solverPressure(n);
					n += 1;
					if (n % 100 == 0)
					{
						Plot_s_emb.Add(s[n - 1]);
						Plot_p_emb.Add(p[n - 1]);
						tempo_plot_emb.Add(t_sim.ToString());
						simLogger(n, t, t_sim);
					}
					t_fim = t_sim;
				}
			}
			embebicao = false;

			Console.WriteLine("TEMPO FINAL EMBEBICAO");
			simLogger(n, t, t_fim);
			Plot_s_emb.Add(s[n - 1]);
			Plot_p_emb.Add(p[n - 1]);
			tempo_plot_emb.Add(t_embebicao.ToString());
		}
		public virtual void solverPressure(int n) {
			double[] a = new double[nbl];
			double[] b = new double[nbl];
			double[] c = new double[nbl];
			double[] d = new double[nbl];
			double[] coefficients;
			for (int i = 0; i < nbl; i++)
			{
				coefficients = preEq(n, i);
				if (i > 0)
				{
					a[i] = coefficients[0];
				}
				if (i < nbl - 1)
				{
					c[i] = coefficients[2];
				}
				b[i] = coefficients[1];
				d[i] = coefficients[3];
			}
			p.Add(thomasSolver(ref a, ref b, ref c, ref d));
		}
		public virtual double solverSaturation1(int n) {
			double[] s_n1 = new double[nbl];
			g = new(nbl);
			g.AddRange(new double[nbl]);
			for (int i = 0; i < nbl; i++)
			{
				satEq2(n, i);
			}
			double t = dt();

			for (int i = 0; i < nbl; ++i)
			{
				s_n1[i] = g[i] * t + s[n - 1][i];
			}
			s.Add(s_n1.ToList());

			return dt();
		}
		public virtual double solverSaturation2(int n) {
			double[] s_n1 = new double[nbl];
			double t = dt();
			g = new(nbl);
			g.AddRange(new double[nbl]);
			double gi;
			for (int i = 0; i < nbl; i++)
			{
				gi = satEq1(n, i);
				s_n1[i] = gi * t + s[n - 1][i];
			}
			s.Add(s_n1.ToList());
			return dt();
		}
		public double[] preEq(int n, int i) {
			double[] c = new double[4];
			//c[0] p_i-1
			//c[1] p_i
			//c[2] p_i+1
			//c[3] lado direito da eq;
			if (i == 0)
			{
				double tsp = upstreamWeighting(n, i, 1);
				double lambda_wPOS = core.lambda_w(tsp, uw);
				double lambda_tPOS = core.lambda_t(tsp, uw, uo);
				c[0] = 0;
				c[1] = -lambda_tPOS;
				c[2] = lambda_tPOS;
				c[3] = ((-u() * h / core.permeability) + lambda_wPOS * (pc[n][i + 1] - pc[n][i]));
			}
			else if (i == nbl - 1)
			{
				double tsp = upstreamWeighting(n, i, 1);
				double tsn = upstreamWeighting(n, i, -1);
				double lambda_wPOS = core.lambda_w(tsp, uw);
				double lambda_wNEG = core.lambda_w(tsn, uw);
				double lambda_tPOS = core.lambda_t(tsp, uw, uo);
				double lambda_tNEG = core.lambda_t(tsn, uw, uo);

				c[0] = lambda_tNEG;
				c[1] = -(2 * lambda_tPOS + lambda_tNEG);
				c[2] = 0;
				c[3] = -2 * lambda_tPOS * Pl - (2 * lambda_wPOS + lambda_wNEG) * pc[n][i] + lambda_wNEG * pc[n][i - 1];
			}
			else
			{
				double tsp = upstreamWeighting(n, i, 1);
				double tsn = upstreamWeighting(n, i, -1);
				double lambda_wPOS = core.lambda_w(tsp, uw);
				double lambda_wNEG = core.lambda_w(tsn, uw);
				double lambda_tPOS = core.lambda_t(tsp, uw, uo);
				double lambda_tNEG = core.lambda_t(tsn, uw, uo);
				c[0] = lambda_tNEG;
				c[1] = -(lambda_tNEG + lambda_tPOS);
				c[2] = lambda_tPOS;
				c[3] = (lambda_wPOS * pc[n][i + 1] - (lambda_wPOS + lambda_wNEG) * pc[n][i] + lambda_wNEG * pc[n][i - 1]);
			}
			return c;
		}
		public double satEq1(int n, int i) {
			double gi = 0.0;
			if (i == 0)
			{
				double Sp = upstreamWeighting(n, i, 1);
				double Ao_POS = core.lambda_o(Sp, uo);
				if (drenagem)
				{
					gi = (-core.permeability / (core.porosity * Math.Pow(h, 2))) * (Ao_POS * (p[n - 1][i + 1] - p[n - 1][i]) + (u() * h / core.permeability));
				}
				else if (embebicao)
				{
					gi = (-core.permeability / (core.porosity * Math.Pow(h, 2))) * (Ao_POS * (p[n - 1][i + 1] - p[n - 1][i]));
				}
			}
			else if (i == nbl - 1)
			{
				double Sp = upstreamWeighting(n, i, 1);
				double Sn = upstreamWeighting(n, i, -1);
				double Ao_POS = core.lambda_o(Sp, uo);
				double Ao_NEG = core.lambda_o(Sn, uo);
				gi = (-4 * core.permeability / (3 * core.porosity * h * h)) * (2 * Ao_POS * (Pl - p[n - 1][i]) - Ao_NEG * (p[n - 1][i] - p[n - 1][i - 1]));

			}
			else
			{
				double Sp = upstreamWeighting(n, i, 1);
				double Sn = upstreamWeighting(n, i, -1);
				double Ao_POS = core.lambda_o(Sp, uo);
				double Ao_NEG = core.lambda_o(Sn, uo);
				gi = (-core.permeability / (core.porosity * Math.Pow(h, 2))) * (Ao_POS * (p[n - 1][i + 1] - p[n - 1][i]) - Ao_NEG * (p[n - 1][i] - p[n - 1][i - 1]));
			}
			g[i] = gi;
			return gi;
		}
		public double satEq2(int n, int i) {
			double gi = 0.0;
			if (i == 0)
			{
				double ts = upstreamWeighting(n, i, 1);
				double fwPOS = core.fw(ts, uw, uo);
				double lamb_oPOS = core.lambda_o(ts, uo);
				double T = fwPOS * lamb_oPOS;
				if (drenagem)
				{
					gi = -1 / core.porosity * (core.permeability / Math.Pow(h, 2) * (T * (pc[n - 1][i + 1] - pc[n - 1][i])) + u() / h * fwPOS);
				}
				else if (embebicao)
				{
					gi = -1 / core.porosity * (core.permeability / Math.Pow(h, 2) * (T * (pc[n - 1][i + 1] - pc[n - 1][i])) + u() / h * (fwPOS - 1));
				}
			}
			else if (i == nbl - 1)
			{
				double tsp = upstreamWeighting(n, i, 1);
				double tsn = upstreamWeighting(n, i, -1);
				double fwPOS = core.fw(tsp, uw, uo);
				double fwNEG = core.fw(tsn, uw, uo);
				double lamb_oPOS = core.lambda_o(tsp, uo);
				double lamb_oNEG = core.lambda_o(tsn, uo);
				double T = fwNEG * lamb_oNEG;
				gi = -1 / core.porosity * ((2 * u() / h) * (fwPOS - fwNEG) + ((4 * core.permeability) / (3 * Math.Pow(h, 2))) * ((2 * fwPOS * lamb_oPOS * (Pl - pc[n - 1][i])) - (T * (pc[n - 1][i] - pc[n - 1][i - 1]))));

			}
			else
			{
				double tsp = upstreamWeighting(n, i, 1);
				double tsn = upstreamWeighting(n, i, -1);
				double fwPOS = core.fw(tsp, uw, uo);
				double fwNEG = core.fw(tsn, uw, uo);
				double lamb_oPOS = core.lambda_o(tsp, uo);
				double lamb_oNEG = core.lambda_o(tsn, uo);
				gi = -1 / core.porosity * (core.permeability / Math.Pow(h, 2) * (lamb_oPOS * fwPOS * pc[n - 1][i + 1] - (lamb_oNEG * fwNEG + lamb_oPOS * fwPOS) * pc[n - 1][i] + lamb_oNEG * fwNEG * pc[n - 1][i - 1]) + u() / h * (fwPOS - fwNEG));
			}
			g[i] = gi;
			return gi;
		}
		public void plotSat(bool multiplot = false) { }
		public void plotPre() { }
		public void grid() {
			h = core.length / nbl;
			for (double k = h; k < nbl; k += h)
			{
				var _h = k - h / 2;
				Plot_h.Append(_h);
			}
			g = new List<double>(new double[nbl]);
		}
		public double upstreamWeighting(int n, int i, int orientation) {
			if (n == 0)
			{
				return s[n][i];
			}
			if (orientation == 1)
			{
				return s[n - 1][i];
			}
			else if (orientation == -1)
			{
				return s[n - 1][i - 1];
			}
			return 0;
		}
		protected void simLogger(int n, double dt, double t_sim)
        {
			Console.WriteLine("Terminou-se o passo de tempo {0} com Dt = {1}\n", n, dt);
			Console.WriteLine("Tempo da simulação total: {0}", t_sim);
		}

		public List<double[]> GetPressureDrainage() { return Plot_p; }
		public List<double[]> GetPressureImbibition() { return Plot_p_emb; }
		public List<List<double>> GetSaturationDrainage() { return Plot_s; }
		public List<List<double>> GetSaturationImbibition() { return Plot_s_emb; }
	}
}
