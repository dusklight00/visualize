"""
geoplot.py
----------

This visualization renders a 3-D plot of the data given the state
trajectory of a simulation, and the path of the property to render.

It generates an HTML file that contains code to render the plot
using Cesium Ion, and the GeoJSON file of data provided to the plot.

The visualization supports different visual representations like color gradients
or size variations to represent the values in the simulation data.

An example of its usage is as follows:

```py
from agent_torch.visualize import GeoPlot

# create a simulation
# ...

# create a visualizer
engine = GeoPlot(config, {
  cesium_token: "...",  # Your Cesium Ion access token
  step_time: 3600,      # Time step in seconds between simulation states
  coordinates = "agents/consumers/coordinates",  # Path to coordinate data in state
  feature = "agents/consumers/money_spent",      # Path to feature data to visualize
  visualization_type: "color",  # Visualization type: "color" or "size"
})

# visualize in the runner-loop
for i in range(0, num_episodes):
  runner.step(num_steps_per_episode)
  engine.render(runner.state_trajectory)
```
"""

import re
import json

import pandas as pd
import numpy as np

from string import Template
from agent_torch.core.helpers import get_by_path

# HTML template with embedded JavaScript for the Cesium visualization
# This template will be populated with simulation data using string.Template
# The $variables will be replaced with actual values during the render process
geoplot_template = """
<!doctype html>
<html lang="en">
	<head>
		<meta charset="UTF-8" />
		<meta
			name="viewport"
			content="width=device-width, initial-scale=1.0"
		/>
		<title>Cesium Time-Series Heatmap Visualization</title>
		<script src="https://cesium.com/downloads/cesiumjs/releases/1.95/Build/Cesium/Cesium.js"></script>
		<link
			href="https://cesium.com/downloads/cesiumjs/releases/1.95/Build/Cesium/Widgets/widgets.css"
			rel="stylesheet"
		/>
		<style>
			#cesiumContainer {
				width: 100%;
				height: 100%;
			}
		</style>
	</head>
	<body>
		<div id="cesiumContainer"></div>
		<script>
			// Your Cesium ion access token here
			Cesium.Ion.defaultAccessToken = '$accessToken'

			// Create the viewer
			const viewer = new Cesium.Viewer('cesiumContainer')

			function interpolateColor(color1, color2, factor) {
				const result = new Cesium.Color()
				result.red = color1.red + factor * (color2.red - color1.red)
				result.green =
					color1.green + factor * (color2.green - color1.green)
				result.blue = color1.blue + factor * (color2.blue - color1.blue)
				result.alpha = '$visualType' == 'size' ? 0.2 :
					color1.alpha + factor * (color2.alpha - color1.alpha)
				return result
			}

			function getColor(value, min, max) {
				const factor = (value - min) / (max - min)
				return interpolateColor(
					Cesium.Color.BLUE,
					Cesium.Color.RED,
					factor
				)
			}

			function getPixelSize(value, min, max) {
				const factor = (value - min) / (max - min)
				return 100 * (1 + factor)
			}

			function processTimeSeriesData(geoJsonData) {
				const timeSeriesMap = new Map()
				let minValue = Infinity
				let maxValue = -Infinity

				geoJsonData.features.forEach((feature) => {
					const id = feature.properties.id
					const time = Cesium.JulianDate.fromIso8601(
						feature.properties.time
					)
					const value = feature.properties.value
					const coordinates = feature.geometry.coordinates

					if (!timeSeriesMap.has(id)) {
						timeSeriesMap.set(id, [])
					}
					timeSeriesMap.get(id).push({ time, value, coordinates })

					minValue = Math.min(minValue, value)
					maxValue = Math.max(maxValue, value)
				})

				return { timeSeriesMap, minValue, maxValue }
			}

			function createTimeSeriesEntities(
				timeSeriesData,
				startTime,
				stopTime
			) {
				const dataSource = new Cesium.CustomDataSource(
					'AgentTorch Simulation'
				)

				for (const [id, timeSeries] of timeSeriesData.timeSeriesMap) {
					const entity = new Cesium.Entity({
						id: id,
						availability: new Cesium.TimeIntervalCollection([
							new Cesium.TimeInterval({
								start: startTime,
								stop: stopTime,
							}),
						]),
						position: new Cesium.SampledPositionProperty(),
						point: {
							pixelSize: '$visualType' == 'size' ? new Cesium.SampledProperty(Number) : 10,
							color: new Cesium.SampledProperty(Cesium.Color),
						},
						properties: {
							value: new Cesium.SampledProperty(Number),
						},
					})

					timeSeries.forEach(({ time, value, coordinates }) => {
						const position = Cesium.Cartesian3.fromDegrees(
							coordinates[0],
							coordinates[1]
						)
						entity.position.addSample(time, position)
						entity.properties.value.addSample(time, value)
						entity.point.color.addSample(
							time,
							getColor(
								value,
								timeSeriesData.minValue,
								timeSeriesData.maxValue
							)
						)

						if ('$visualType' == 'size') {
						  entity.point.pixelSize.addSample(
  							time,
  							getPixelSize(
  								value,
  								timeSeriesData.minValue,
  								timeSeriesData.maxValue
  							)
  						)
						}
					})

					dataSource.entities.add(entity)
				}

				return dataSource
			}

			// Example time-series GeoJSON data
			const geoJsons = $data

			const start = Cesium.JulianDate.fromIso8601('$startTime')
			const stop = Cesium.JulianDate.fromIso8601('$stopTime')

			viewer.clock.startTime = start.clone()
			viewer.clock.stopTime = stop.clone()
			viewer.clock.currentTime = start.clone()
			viewer.clock.clockRange = Cesium.ClockRange.LOOP_STOP
			viewer.clock.multiplier = 3600 // 1 hour per second

			viewer.timeline.zoomTo(start, stop)

			for (const geoJsonData of geoJsons) {
				const timeSeriesData = processTimeSeriesData(geoJsonData)
				const dataSource = createTimeSeriesEntities(
					timeSeriesData,
					start,
					stop
				)
				viewer.dataSources.add(dataSource)
				viewer.zoomTo(dataSource)
			}
		</script>
	</body>
</html>
"""


def read_var(state, var):
    """
    Access a nested property in the state object using path notation.
    
    Args:
        state (dict): The state object from which to extract data
        var (str): Path to the variable in format "path/to/property"
        
    Returns:
        The value at the specified path in the state dictionary
    """
    return get_by_path(state, re.split("/", var))


class GeoPlot:
    """
    GeoPlot visualizer for AgentTorch simulations.
    
    This class creates interactive 3D map visualizations using Cesium Ion to display 
    time-series data from simulation state trajectories. It generates both a GeoJSON 
    data file and an HTML file with embedded JavaScript for the interactive visualization.
    
    The visualization can represent values using either color gradients (from blue to red)
    or varying dot sizes, based on the visualization_type parameter.
    """
    
    def __init__(self, config, options):
        """
        Initialize the GeoPlot visualizer.
        
        Args:
            config (dict): Configuration dictionary with simulation metadata
            options (dict): Visualization options containing:
                - cesium_token (str): Cesium Ion API access token
                - step_time (int): Time in seconds between simulation steps
                - coordinates (str): Path to coordinates data in state (format: "path/to/coordinates")
                - feature (str): Path to feature data to visualize (format: "path/to/feature")
                - visualization_type (str): Type of visualization - "color" or "size"
        """
        self.config = config
        (
            self.cesium_token,
            self.step_time,
            self.entity_position,
            self.entity_property,
            self.visualization_type,
        ) = (
            options["cesium_token"],
            options["step_time"],
            options["coordinates"],
            options["feature"],
            options["visualization_type"],
        )

    def render(self, state_trajectory):
        """
        Generate visualization from simulation state trajectory.
        
        This method processes the simulation state trajectory to extract coordinates and
        feature values, then generates both a GeoJSON data file and an HTML visualization file.
        
        Args:
            state_trajectory (list): List of state sequences from simulation episodes
                                     Format: [episode_states_1, episode_states_2, ...]
                                     where each episode_states is a list of state objects
        
        Returns:
            None: Files are written to disk with names based on simulation metadata
        """
        # Initialize empty containers for coordinates and values
        coords, values = [], []
        name = self.config["simulation_metadata"]["name"]
        geodata_path, geoplot_path = f"{name}.geojson", f"{name}.html"

        # Extract coordinates and feature values from each episode's final state
        for i in range(0, len(state_trajectory) - 1):
            final_state = state_trajectory[i][-1]  # Get final state of each episode

            # Extract and convert coordinate and feature data
            coords = np.array(read_var(final_state, self.entity_position)).tolist()
            values.append(
                np.array(read_var(final_state, self.entity_property)).flatten().tolist()
            )

        # Generate timestamps for visualization timeline
        start_time = pd.Timestamp.utcnow()
        timestamps = [
            start_time + pd.Timedelta(seconds=i * self.step_time)
            for i in range(
                self.config["simulation_metadata"]["num_episodes"]
                * self.config["simulation_metadata"]["num_steps_per_episode"]
            )
        ]

        # Create GeoJSON features for each coordinate and timestamp
        geojsons = []
        for i, coord in enumerate(coords):
            features = []
            for time, value_list in zip(timestamps, values):
                # Create GeoJSON feature with point geometry and properties
                features.append(
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "Point",
                            # Note: GeoJSON uses [longitude, latitude] order
                            "coordinates": [coord[1], coord[0]],
                        },
                        "properties": {
                            "value": value_list[i],  # Value to visualize
                            "time": time.isoformat(),  # ISO format timestamp
                        },
                    }
                )
            # Create a GeoJSON feature collection for each coordinate
            geojsons.append({"type": "FeatureCollection", "features": features})

        # Write the GeoJSON data to file
        with open(geodata_path, "w", encoding="utf-8") as f:
            json.dump(geojsons, f, ensure_ascii=False, indent=2)

        # Create HTML visualization from template
        tmpl = Template(geoplot_template)
        with open(geoplot_path, "w", encoding="utf-8") as f:
            f.write(
                tmpl.substitute(
                    {
                        "accessToken": self.cesium_token,
                        "startTime": timestamps[0].isoformat(),
                        "stopTime": timestamps[-1].isoformat(),
                        "data": json.dumps(geojsons),
                        "visualType": self.visualization_type,
                    }
                )
            )
