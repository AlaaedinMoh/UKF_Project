/* \author Aaron Brown */
// Create simple 3d highway enviroment using PCL
// for exploring self-driving car sensors

//#include "render/render.h"
#include "highway.h"

int main(int argc, const char* argv[])
{
	bool useLaser = true;
	bool useRdr = true;
	std::vector<bool> showCars = {true, true, true};
	if(argc > 1)
	{
		string argStr = argv[1];
	    useLaser = argStr.compare("l")==0 || argStr.compare("lr")==0;	
		useRdr = argStr.compare("r")==0 || argStr.compare("lr")==0;
		cout<<"Args 2 = "<<argStr<<endl; 
		if(argc > 2)
		{
			argStr = argv[2];
			showCars[0] = argStr.compare("t")==0;
			if(argc > 3)
			{
				argStr = argv[3];
				showCars[1] = argStr.compare("t")==0;
				if(argc > 4)
				{
					argStr = argv[4];
					showCars[2] = argStr.compare("t")==0;
				}
			}
		}
	}
	cout<<"Visualize Lidar = "<<std::to_string(useLaser)<<endl;
	cout<<"Visualize Radar = "<<std::to_string(useRdr)<<endl;
	cout<<"Show first car = "<<std::to_string(showCars[0])<<endl;
	cout<<"Show second car = "<<std::to_string(showCars[1])<<endl;
	cout<<"Show second car = "<<std::to_string(showCars[2])<<endl;
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);

	// set camera position and angle
	viewer->initCameraParameters();
	float x_pos = 0;
	viewer->setCameraPosition ( x_pos-26, 0, 15.0, x_pos+25, 0, 0, 0, 0, 1);

	Highway highway(viewer);
	highway.visualize_lidar = useLaser;
	highway.visualize_radar = useRdr;
	highway.trackCars = showCars;

	//initHighway(viewer);

	int frame_per_sec = 30;
	int sec_interval = 10;
	int frame_count = 0;
	int time_us = 0;

	double egoVelocity = 25;

	while (frame_count < (frame_per_sec*sec_interval))
	{
		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

		//stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		highway.stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		viewer->spinOnce(1000/frame_per_sec);
		frame_count++;
		time_us = 1000000*frame_count/frame_per_sec;
		
	}

}