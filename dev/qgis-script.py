from qgis.core import (QgsRasterLayer, QgsPointXY, QgsProject,
                       QgsVectorLayer, QgsCoordinateTransform,
                       QgsCoordinateReferenceSystem, QgsField,
                       QgsPointXY)
from PyQt5.QtCore import QVariant


def cartesian_distance(row1, col1, row2, col2, pixel_width, pixel_height):
    # Calculate the Cartesian distance between two cells in a raster grid
    dx = (col1 - col2) * pixel_width
    dy = (row1 - row2) * pixel_height
    return (dx**2 + dy**2)**0.5


def find_nearest_cell_with_data(raster_layer, point, max_distance):
    # Function to check if the given cell is valid and has data
    def is_valid_and_has_data(row, col):
        if 0 <= row < raster_layer.height() and 0 <= col < raster_layer.width():
            value = raster_layer.dataProvider().sample(QgsPointXY(col, row), 1)[0]
            return value is not None and value != raster_layer.dataProvider().sourceNoDataValue(1)
        return False

    # Convert the input point to raster cell coordinates
    start_row, start_col = get_raster_cell(raster_layer, point)

    # Start with the initial cell
    if is_valid_and_has_data(start_row, start_col):
        return start_row, start_col

    pixel_width = raster_layer.rasterUnitsPerPixelX()
    pixel_height = raster_layer.rasterUnitsPerPixelY()

    # Spiral search pattern
    distance_threshold = max_distance * max(pixel_width, pixel_height)
    for d_row in range(-max_distance, max_distance + 1):
        for d_col in range(-max_distance, max_distance + 1):
            row, col = start_row + d_row, start_col + d_col
            if is_valid_and_has_data(row, col):
                # Calculate Cartesian distance
                distance = cartesian_distance(start_row, start_col, row, col, pixel_width, pixel_height)
                if distance <= distance_threshold:
                    return row, col

    # No valid cell found within max_distance
    return None, None


def get_raster_cell(raster_layer, point):
    # Get the extent of the raster layer
    extent = raster_layer.renderer().extent()

    # Check if the point is within the extent of the raster layer
    if extent.contains(point):
        # Calculate the pixel coordinates (column and row) based on the extent and point
        x = (point.x() - extent.xMinimum()) / raster_layer.rasterUnitsPerPixelX()
        y = (extent.yMaximum() - point.y()) / raster_layer.rasterUnitsPerPixelY()

        # Convert to integers to get the row and column
        col = int(x)
        row = int(y)

        return row, col

    return None, None  # Point is outside the raster extent


def get_neighbor_values(raster_layer, row, col, max_distance):
    neighbor_values = []

    # Get the number of rows and columns in the raster layer
    num_rows = raster_layer.height()
    num_cols = raster_layer.width()

    # Get the band count
    band_count = raster_layer.bandCount()

    # Determine the range of neighboring rows and columns
    for i in range(-max_distance, max_distance + 1):
        for j in range(-max_distance, max_distance + 1):
            neighbor_row = row + i
            neighbor_col = col + j

            # Check if the neighbor cell is within the raster bounds
            if (0 <= neighbor_row < num_rows) and (0 <= neighbor_col < num_cols):
                for band_index in range(1, band_count + 1):
                    # Convert the row and column back to QgsPointXY
                    neighbor_point = QgsPointXY(
                        neighbor_col * raster_layer.rasterUnitsPerPixelX() + raster_layer.extent().xMinimum(),
                        raster_layer.extent().yMaximum() - neighbor_row * raster_layer.rasterUnitsPerPixelY()
                    )

                    # Get the value of the neighbor cell for each band
                    value, success = raster_layer.dataProvider().sample(neighbor_point, band_index)

                    if value is not None and value != raster_layer.dataProvider().sourceNoDataValue(band_index):
                        neighbor_values.append(value)

    return neighbor_values


def sample_mean_neighbors(raster_layer, point, max_distance):
    # Convert the point to its corresponding raster cell (row, column)
    row, col = get_raster_cell(raster_layer, point)

    # Get the values of neighboring cells
    neighbor_values = get_neighbor_values(raster_layer, row, col, max_distance)

    # Filter out invalid data and calculate the mean
    valid_values = [value for value in neighbor_values if value is not None]
    if valid_values:
        mean_value = sum(valid_values) / len(valid_values)
        return round(mean_value, 3)  # Round to 3 decimal places
    else:
        return "N/A"


def sample_raster_at_point(raster_layer, point):
    # Get raster data provider
    provider = raster_layer.dataProvider()

    # Sample the raster at the exact point
    value, success = provider.sample(point, 1)

    # Check if the sampling was successful
    if success and not value is None:
        return round(value, 3)  # Round to 3 decimal places
    else:
        return "N/A"


def sample_nearest_cell(raster_layer, point, max_distance):    
    # Get the nearest cell with data
    nearest_row, nearest_col = find_nearest_cell_with_data(raster_layer, point, max_distance)

    if nearest_row is not None:
        # Convert the nearest cell coordinates to QgsPointXY
        nearest_point = QgsPointXY(
            nearest_col * raster_layer.rasterUnitsPerPixelX() + raster_layer.extent().xMinimum(),
            raster_layer.extent().yMaximum() - nearest_row * raster_layer.rasterUnitsPerPixelY()
        )

        # Sample the value at the nearest cell
        value = sample_raster_at_point(raster_layer, nearest_point)

        if value is not None:
            return round(value, 3)  # Round to 3 decimal places

    return "N/A"


def perform_raster_sampling(shapefile_path, raster_file_path, attributive_field_name, method, max_distance=1):
    # Ensure the attributive field name is not longer than 8 characters
    if len(attributive_field_name) > 8:
        attributive_field_name = attributive_field_name[:8]

    # Load the shapefile as a vector layer
    vector_layer = QgsVectorLayer(shapefile_path, "Sample Points", "ogr")

    # Load the raster file as a raster layer
    raster_layer = QgsRasterLayer(raster_file_path, "Raster Layer")

    # Ensure both layers are loaded
    if not vector_layer.isValid():
        print("Vector layer failed to load!")
        return

    if not raster_layer.isValid():
        print("Raster layer failed to load!")
        return

    # Check if the layers are already added to the project
    if vector_layer not in QgsProject.instance().mapLayers().values():
        # Add the vector layer to the project for visualization
        QgsProject.instance().addMapLayer(vector_layer)

    if raster_layer not in QgsProject.instance().mapLayers().values():
        # Add the raster layer to the project for visualization
        QgsProject.instance().addMapLayer(raster_layer)

    # Check if the Coordinate Reference Systems (CRS) of the vector and raster layers are the same.
    crs_src = vector_layer.crs()
    crs_dest = raster_layer.crs()

    if not crs_src == crs_dest:
        print("Error: Coordinate systems of the vector and raster layers do not match.")
        return

    # Coordinate transformation setup
    transform = QgsCoordinateTransform(crs_src, crs_dest, QgsProject.instance())

    # Check if the attributive field already exists
    if attributive_field_name not in [field.name() for field in vector_layer.fields()]:
        # Field does not exist, so let's create it
        field = QgsField(attributive_field_name, QVariant.Double)
        vector_layer.dataProvider().addAttributes([field])
        vector_layer.updateFields()

    # Index of the new field
    field_index = vector_layer.fields().indexOf(attributive_field_name)

    # Start editing the vector layer
    vector_layer.startEditing()

    # Loop through each feature in the shapefile
    for feature in vector_layer.getFeatures():
        geom = feature.geometry()

        # Check if the geometry is a point
        if geom.type() == QgsWkbTypes.PointGeometry:
            point = geom.asPoint()

            # Transform point to raster CRS
            point_transformed = transform.transform(point)

            # Sample the raster based on the method
            if method == 'exact':
                value = sample_raster_at_point(raster_layer, point_transformed)
            elif method == 'nearest':
                value = sample_nearest_cell(raster_layer, point_transformed, max_distance)
            elif method == 'mean':
                value = sample_mean_neighbors(raster_layer, point_transformed, max_distance)

            # Store the sampled value in the attribute table
            vector_layer.changeAttributeValue(feature.id(), field_index, value)

        else:
            print(f"Skipping non-point geometry: {geom.type()}")

    # Save changes to the layer
    vector_layer.commitChanges()

    # Save the changes to an existing shapefile
    QgsVectorFileWriter.writeAsVectorFormat(vector_layer, '/Users/abegovic4/Desktop/raster_value_sampler/dev/Output_Boulders.shp', 'UTF-8', vector_layer.crs(), 'ESRI Shapefile')

    # Get and print the column names
    field_names = [field.name() for field in vector_layer.fields()]
    print("Column Names:", field_names)

    # Iterate through the features to print attributes
    print("\nFeature Attributes:")
    for feature in vector_layer.getFeatures():
        attrs = feature.attributes()
        print(dict(zip(field_names, attrs)))


shapefile_path = '/Users/abegovic4/Desktop/raster_value_sampler/dev/Output_Boulders.shp'
raster_file_path = '/Users/abegovic4/Desktop/raster_value_sampler/dev/Test_Encoded_Depths_File.tif'

# Example call of the main function
perform_raster_sampling(
    shapefile_path=shapefile_path,
    raster_file_path=raster_file_path,
    attributive_field_name='RasterVal',
    method='nearest',
    max_distance=3
)
