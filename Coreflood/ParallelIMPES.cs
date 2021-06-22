using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Coreflood
{
	public class ParallelIMPES : IMPES
	{
		bool isPR1;
		int[] wkload;
		int num_of_threads;

		List<int> dre;
		List<int> emb;
		List<double> s_n1;
		List<double> pc0;

		static Mutex m = new();
		Thread[] threads;
		ManualResetEvent[] handles;

		public ParallelIMPES(RCore ncore, IMPESArgs ia, bool _isPR1): base(ncore, ia){
			num_of_threads = Environment.ProcessorCount;
			threads = new Thread[num_of_threads];
			wkload = workload(nbl, num_of_threads);
			isPR1 = _isPR1;
			dre = new();
			emb = new();
		}

		public override void solverManifold(){
			double[] s0 = new double[nbl];
			double[] pc0 = new double[nbl];

			Array.Fill(s0, si);
			Array.Fill(pc0, core.Pct);

			s.Add(s0.ToList());
			Plot_s.Add(s0.ToList());
			tempo_plot.Add("0");
			pc.Add(pc0.ToList());

			int n = 1;
			double t = 0;

			solveDrainage(ref n, ref t);
			solveImbibition(ref n, ref t);

			Console.WriteLine("Resolvendo o campo da pressao...\n");
			parallelPressureSolver(workload(s.Count, num_of_threads));
			Console.WriteLine("Simulacao completa!\n");

			foreach (var pos in dre)
			{
				Plot_s.Add(s[pos]);
				Plot_p.Add(p[pos]);
			}

			foreach (var pos in emb)
			{
				Plot_s_emb.Add(s[pos]);
				Plot_p_emb.Add(p[pos]);
			}
		}
		private void solveDrainage(ref int n, ref double t) {
			drenagem = true;
			for (double t_sim = 0; t_sim < t_drenagem; t_sim += t)
			{
				t = solverSaturation1(n);
				n += 1;
				if (n % 2000 == 0)
				{
					dre.Add(n - 1);
				}
				simLogger(n, t, t_sim);
			}
			dre.Add(n - 1);
			drenagem = false;
		}
		private void solveImbibition(ref int n, ref double t) {
			embebicao = true;
			for (double t_sim = t_drenagem; t_sim < t_embebicao; t_sim += t)
			{
				t = solverSaturation1(n);
				n += 1;
				if (n % 100 == 0)
				{
					emb.Add(n - 1);
				}
				simLogger(n, t, t_sim);
			}
			emb.Add(n - 1);
			embebicao = false;
		}
		public void solverPressure(int start, int end) {
			for (int n = start; n < end; n++)
			{
				double[] a = new double[nbl];
				double[] b = new double[nbl]; 
				double[] c = new double[nbl]; 
				double[] d = new double[nbl];
				double[] coefficients = new double[4];
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
				p[n] = thomasSolver(ref a, ref b, ref c, ref d);
		}
	}
		public override double solverSaturation1(int n) {
			g = new(nbl);
			s_n1 = new(nbl);
			pc0 = new(nbl);

			g.AddRange(new double[nbl]);
			s_n1.AddRange(new double[nbl]);
			pc0.AddRange(new double[nbl]);

			if (isPR1)
			{
				for (int i = 0; i < nbl; i++)
				{
					satEq2(n, i);
				}
			}
			else
			{
				handles = new ManualResetEvent[num_of_threads];
				for (int i = 0; i < num_of_threads; i++)
				{
					handles[i] = new ManualResetEvent(false);
					ThreadPool.QueueUserWorkItem(ThreadPoolCallback, new Tuple<int, int>(i, n));
				}
				WaitHandle.WaitAll(handles);
			}
			double t = dt();

			for (int i = 0; i < nbl; ++i)
			{
				s_n1[i] = g[i] * t + s[n - 1][i];
				pc0[i] = core.pc(s_n1[i]);
			}
			s.Add(s_n1);
			pc.Add(pc0);
			return t;
		}
		private void ThreadPoolCallback(object state)
        {
			var arg_pack = (Tuple<int, int>)state;
			int i = arg_pack.Item1;
			int n = arg_pack.Item2;

			var s1 = wkload[i];
			var s2 = wkload[i+1];
			satIncrementParallel(s1, s2, n);

			handles[i].Set();
		}
		public void parallelPressureSolver(int[] wkload) {
			p = new double[s.Count][].ToList();
			
			for (int i = 0; i < num_of_threads; ++i)
			{
				var s1 = wkload[i];
				var s2 = wkload[i + 1];
				threads[i] = new Thread(() => solverPressure(s1, s2));
				threads[i].Start();
			}
			foreach (var thread in threads)
			{
				thread.Join();
			}
		}
		public void satIncrementParallel(int start, int end, int n) {
			for (int i = start; i < end; i++)
			{
				satEq2(n, i);
			}
		}
		public static int[] workload(int tasks, int num_of_workers) {
			int remainder = tasks % num_of_workers;
			int equal_tasks = (tasks - remainder) / num_of_workers;

			List<int> tasks_for_each_thread = new int[num_of_workers].ToList();
			int last = 0;
			for (int j = 0; j < num_of_workers; j++)
			{
				tasks_for_each_thread[j] = last + equal_tasks;
				last += equal_tasks;
			}
			for (int i = 0; i < remainder; i++)
			{
				tasks_for_each_thread[i] += 1;
				for (int k = i + 1; k < tasks_for_each_thread.Count; k++)
				{
					tasks_for_each_thread[k] += 1;
				}
			}
			tasks_for_each_thread.Insert(0, 0);
			return tasks_for_each_thread.ToArray();
		}
	}
}
