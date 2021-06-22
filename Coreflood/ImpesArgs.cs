using System;
using System.ComponentModel.DataAnnotations;

namespace Coreflood
{
	public class IMPESArgs
	{
		[Required]
		public double atmospheric_pressure {get; set;}
		[Required]
		public double injection_flow {get; set;}
		[Required]
		public double oil_viscosity {get; set;}
		[Required]
		public double water_viscosity{get; set;}
		[Required]
		public double maximum_saturation_step {get; set;}
		[Required]
		public double pore_volume_injected_in_drainage {get; set;}
		[Required]
		public double pore_volume_injected_in_imbibition {get; set;}
		[Required]
		[Range(0,1, ErrorMessage = "Absolute saturation values are between 0 and 1.")]
		public double initial_saturation {get; set;}
		[Required]
		[Range(1, 1000,ErrorMessage = "Insert a number between 0 and 1000.")]
		public int number_of_blocks {get; set;}
	};
}
