using System;
using System.ComponentModel.DataAnnotations;

namespace Coreflood
{
	public class CoreArgs
	{
		[Required]
		[Range(0.0, 1.0, ErrorMessage = "Porosity Values must be between 0 and 1")]
		public double porosity { get; set; }
		[Required]
		public double diameter {get; set;} 
		[Required]
		public double length {get; set;}
		[Required]
		public double permeability {get; set;}
		[Required]
		public double exp_krw {get; set;}
		[Required]
		public double exp_kro {get; set;}
		[Required]
		public double exp_pc {get; set;}
		[Required]
		public double krwMax {get; set;}
		[Required]
		public double kroMax {get; set;}
		[Required]
		public double pcMax {get; set;}
		[Required]
		[Range(0, 1, ErrorMessage = "Saturation Values must be between 0 and 1")]
		public double sor {get; set;}
		[Required]
		[Range(0, 1, ErrorMessage = "Saturation Values must be between 0 and 1")]
		public double swi {get; set;}
		[Required]
		public double pct {get; set;}
	};
}
