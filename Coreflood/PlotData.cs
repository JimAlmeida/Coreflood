using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coreflood
{
    public class PlotData
    {
        public List<double[]> drainage_pressure { get; set; }
        public List<double[]> imbibition_pressure { get; set; }
        public List<List<double>> drainage_saturation { get; set; }
        public List<List<double>> imbibition_saturation { get; set; }
        public int blocks { get; set;} //each value in the x-axis corresponds to one block as if the chart area was the entire core plug segmented by each one of the blocks.
        public List<string> time_labels { get; set; }

    }
}
