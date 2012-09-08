# ---------------------------------------------------------------------------
# opendatamap.py
# ---------------------------------------------------------------------------
# Generate a heat map from spatial data in a CKAN instance.
# Currently only works with shapefiles.
#----------------------------------------------------------------------------
# This script completes the following:
# 1) Accesses a CKAN instance and downloads all datasets that contain 
#    shapefile resources
# 2) Extracts shapefiles and re-project to web mercator
# 3) Rasterize shapefiles using GDAL and PIL
# 4) Generate a composite image file from all rasterized datasets
# 5) Colorize the composite image file
# 6) The output image file can be used directly in gdal2tiles to 
#    create map tiles for use on a web map
# ---------------------------------------------------------------------------
#----------------------------------------------------------------------------
import os, sys, Image, ImageFilter, urllib2, zipfile, shutil, numpy, math
import ckanclient
from osgeo import gdal, ogr, osr

#TODO: Implement argparse for passing arguments

# Globals
datasets_folder = "datasets"
download_folder = "download"

ckan_host = "http://data.opencolorado.org/api/2"

# Bounding box in web mercator coordinates
x_min = -12140532.1637
x_max = -11359138.5791
y_min = 4438050.84302
y_max = 5012849.66619

# Pixel size in web mercator meters
resolution = float(76.43702)

force_rasterize = False

def main():
    
    print "Running opendatamap.py"
    
    initialize()
    
    process_ckan_datasets()
    
    generate_composite_image()
    
    colorize_composite_image()
    
    generate_world_file()
    
    print "Completed opendatamap.py"

def initialize():
    
    global datasets_folder
    
    # Create the datasets download and analysis folder if it does not exist
    if not os.path.exists(datasets_folder):
        os.makedirs(datasets_folder)
      
def process_ckan_datasets():
    
    global ckan_host
    
    # Initialize the CKAN client  
    ckan_client = ckanclient.CkanClient(base_location=ckan_host)
    
    package_id_list = ckan_client.package_register_get()
    
    index = 0;
    for package_id in package_id_list:
            
        # Get the package details
        package = ckan_client.package_entity_get(package_id)
        
        # Get the package name (slug)
        package_name = package['name']
        
        print "------------------------------"
        print "Processing dataset " + str(index) + " of " + str(len(package_id_list)) + ": " + package_name
        print "Created: " + package['metadata_created'] + ", modified: " + package['metadata_modified']
        
        shapefile_found = False
        resources = package['resources']
        for resource in resources:
            
            # Look for a shapefile resource
            if (resource['mimetype'] and 'shp' in resource['mimetype'].lower()) or \
               (resource['mimetype_inner'] and 'shp' in resource['mimetype_inner'].lower()) or \
               (resource['format'] and 'shp' in resource['format'].lower()) or \
               (resource['format'] and 'shapefile' in resource['format'].lower()) or \
               (resource['name'] and 'shp' in resource['name'].lower()) or \
               (resource['name'] and 'shapefile' in resource['name'].lower()) or \
               (resource['description'] and 'shp' in resource['description'].lower()) or \
               (resource['description'] and 'shapefile' in resource['description'].lower()):
                
                shapefile_found = True
                process_package_resource(package, resource)
        
        if shapefile_found == False:
            print "No shapefile found."
                                            
        index = index + 1    

def initialize_dataset_folder(package_name):
    global datasets_folder, download_folder
    
    dataset_folder = os.path.join(datasets_folder,package_name)
    dataset_download_folder = os.path.join(dataset_folder, download_folder)
    
    if not os.path.exists(dataset_folder):
        os.makedirs(dataset_folder)
        
    if not os.path.exists(dataset_download_folder):
        os.makedirs(dataset_download_folder)
    
    return dataset_folder

def process_package_resource(package, resource):
    
    global datasets_folder, force_rasterize,download_folder
            
    # Get the package name
    package_name = package["name"]
    
    # Get the resource URL
    url = resource["url"]
    
    # Get the dataset folder
    dataset_folder = os.path.join(datasets_folder,package_name)
                        
    # Check if the map image exists for this dataset before downloading again
    # TODO: Compare timestamp of dataset and image and compare image size (may want to regen at different size)
    if force_rasterize or not os.path.exists(os.path.join(dataset_folder,"map.png")):
    
        # If the projected shapefile exists don't download it again
        shapefile_projected = os.path.join(dataset_folder,"download","projected","projected.shp")
        if not os.path.exists(shapefile_projected):
            # Download the shapefile
            shapefile = download_shapefile(package_name, url)

            # Reproject the shapefile
            if (shapefile != None):
                projected_shapefile = reproject_shapefile(package_name, shapefile)
        else:
            shapefile = shapefile_projected
        
        # Rasterize the projected shapefile
        if os.path.exists(shapefile_projected):
            rasterize_shapefile(package_name, shapefile_projected)
            
            # Delete the download folder
            # Delete the source shapefile
            try:
                shutil.rmtree(os.path.join(package_name,download_folder))
            except:
                print "Unable to delete the downloaded shapefile.  Skipping..."
                
        else:
            print "Unable to find shapefile in downloaded resource: " + url
    else:
        print "Map image exists, skipping dataset.."

    
def download_shapefile(package_name,url):
    global download_folder
    
    shapefile = None
    
    # Initialize the folder to hold the dataset
    dataset_folder = initialize_dataset_folder(package_name)
    
    dataset_download_folder = os.path.join(dataset_folder, download_folder)
    dataset_download_folder_source = os.path.join(dataset_download_folder,"source")
    
    dataset_download_file = os.path.join(dataset_download_folder, "shapefile.zip")
    
    try:
        print "Downloading.."
        request = urllib2.urlopen(url)
        with open(dataset_download_file, 'wb') as fp:
            shutil.copyfileobj(request, fp)
        
        # Unzip the file
        print "Unzipping.."
        zip = zipfile.ZipFile(dataset_download_file)
        zip.extractall(dataset_download_folder_source)
        zip.close()
    
        # Delete the zip file
        os.remove(dataset_download_file)
        
        # Find the shapefile (if this zip actually contains a shapefile)
        for (path, dirs, files) in os.walk(dataset_download_folder_source):
            for file in files:
                if file.endswith(".shp"):
                    shapefile = os.path.join(path,file)
                    break
    except:
        print "Error downloading file: "+ url
        
    return shapefile

def reproject_shapefile(package_name, shapefile):
    
    global download_folder
    
    print "Reprojecting shapefile..."
    
    dataset_folder = initialize_dataset_folder(package_name)
    dataset_download_folder = os.path.join(dataset_folder, download_folder)
    
    source_folder = os.path.join(dataset_download_folder,"source")
    projected_folder = os.path.join(dataset_download_folder,"projected")
    projected_shapefile = os.path.join(projected_folder, "projected.shp")
    projected_shapefile_prj = os.path.join(projected_folder, "projected.prj")
    
    if not os.path.exists(projected_folder):
        os.makedirs(projected_folder)
        
    driver = ogr.GetDriverByName('ESRI Shapefile')
     
    src_shapefile = ogr.Open(shapefile)
    
    src_layer = src_shapefile.GetLayer()
    
    # Get the input SpatialReference
    src_sr = src_layer.GetSpatialRef()

    # create the output SpatialReference
    dest_sr = osr.SpatialReference()
    dest_sr.ImportFromEPSG(3857)
    
    # create the CoordinateTransformation
    transformation = osr.CoordinateTransformation(src_sr, dest_sr)
    
    # create a new data source and layer
    if os.path.exists(projected_shapefile):
        driver.DeleteDataSource(projected_shapefile)
      
    dest_shapefile = driver.CreateDataSource(projected_shapefile)
    
    if dest_shapefile is None:
        print 'Could not create file'
        sys.exit(1)
        
    dest_layer = dest_shapefile.CreateLayer('output', geom_type=src_layer.GetGeomType())
    
    # get the layer definition for the output shapefile
    dest_layer_defn = dest_layer.GetLayerDefn()
    
    # Get the first source feature
    src_feature = src_layer.GetNextFeature()
    while src_feature:
    
        # get the input geometry
        src_geom = src_feature.GetGeometryRef()
        
        if (src_geom != None): 
            # reproject the geometry
            src_geom.Transform(transformation)
            
            # create a new feature
            dest_feature = ogr.Feature(dest_layer_defn)
            
            # set the geometry and attribute
            dest_feature.SetGeometry(src_geom)
            
            # add the feature to the shapefile
            dest_layer.CreateFeature(dest_feature)
            
            # destroy the features and get the next input feature
            dest_feature.Destroy
            src_feature.Destroy
            
            src_feature = src_layer.GetNextFeature()
        else:
            projected_shapefile = None
            print "Unable to load source geometry"
            break
    
    # close the shapefiles
    dest_shapefile.Destroy()
    src_shapefile.Destroy()
    
    # Delete the source shapefile
    try:
        shutil.rmtree(source_folder)
    except:
        print "Unable to delete the source shapefile (" + source_folder + ").  Skipping..."
        
    # create the *.prj file
    if projected_shapefile != None:
        file = open(projected_shapefile_prj, 'w')
        file.write('PROJCS["WGS_1984_Web_Mercator_Auxiliary_Sphere",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Mercator_Auxiliary_Sphere"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],PARAMETER["Auxiliary_Sphere_Type",0.0],UNIT["Meter",1.0]]')
        file.close()  
    
    return projected_shapefile

def rasterize_shapefile(package_name, shapefile):
    """Rasterizes a shapefile to be used in building the composite
    heat map
    """
    global resolution, x_min, x_max, y_min, y_max
    
    print "Rasterizing... "
                 
    projected_shapefile = ogr.Open(shapefile)
    
    layer = projected_shapefile.GetLayer()
        
    # Get the spatial reference
    srs = layer.GetSpatialRef()
            
    # Create the destination data source
    x_res = int((x_max - x_min) / resolution)
    y_res = int((y_max - y_min) / resolution)
        
    output = os.path.join(datasets_folder, package_name, 'map.tif')
    output_png = os.path.join(datasets_folder, package_name, 'map.png')

    target_ds = gdal.GetDriverByName('GTiff').Create(output, x_res, y_res, 3, gdal.GDT_Byte)
    
    target_ds.SetGeoTransform((
            x_min, resolution, 0,
            y_max, 0, -resolution
      ))
    
    # Set the project of the target data source
    if srs:
        # Make the target raster have the same projection as the source
        target_ds.SetProjection(srs.ExportToWkt())
        
    # Set the burn value based on the geometry type
    geom_type=layer.GetGeomType()
    if (geom_type == ogr.wkbPolygon):
        
        # Load the shapefile into memory for modification
        projected_shapefile_memory = ogr.GetDriverByName("Memory").CopyDataSource(projected_shapefile, "")
        layer_memory = projected_shapefile_memory.GetLayer()
        
        # TODO: For polygons, use the polygon size to determine the burn value
        # The basic idea is that large polygons shouldn't have a big influence on 
        # open data availability as there isn't a high data density when the polygons are 
        # large (ex. statewide districts)
        
        # Create a field in the source layer to hold the burn value
        burn_field = "burn_value"
        field_def = ogr.FieldDefn(burn_field, ogr.OFTReal)
        layer_memory.CreateField(field_def)
        layer_memory_def = layer_memory.GetLayerDefn()
        
        field_index = layer_memory_def.GetFieldIndex(burn_field)
        
        # Generate random values for the color field (it's here that the value
        # of the attribute should be used, but you get the idea)
        for feature in layer_memory:
            geometry = feature.GetGeometryRef()
            
            polygon_area = geometry.GetArea()
            
            burn_value = compute_polygon_burn_value(polygon_area)
            
            feature.SetField(field_index, burn_value)
            layer_memory.SetFeature(feature)
            
            geometry = None
            feature = None
            
        gdal.RasterizeLayer(target_ds, (1, 2, 3), layer_memory, None, options=["ATTRIBUTE=%s" % burn_field])
        
        layer_memory = None
        projected_shapefile_memory = None
        
    else:
        gdal.RasterizeLayer(target_ds, (1, 2, 3), layer, burn_values=(255, 255, 255), options=[])
    
    target_ds.FlushCache()
    target_ds = None
    
    # Open the tif image
    im = Image.open(output)
    
    # Convert to 8 bit
    im = im.convert("L")
    
    # If the type is not polygon, apply some smoothing
    #if (geom_type != ogr.wkbPolygon):
    #    im = im.filter(ImageFilter.SMOOTH)
    
    # Save to a PNG image
    im.save(output_png)
    
    print "Created map image for dataset: " + package_name
    
    # Delete the TIF
    os.remove(output)

def generate_composite_image():
    
    global datasets_folder
    
    image_count = 0
    base_image_array = None
    
    print "------------------------------"
    print "Generating composite image"
    

    # Find the image maps for each dataset
    for folder in os.listdir(datasets_folder):
        for file in os.listdir(os.path.join(datasets_folder,folder)):
            if file == "map.png":
             
                print "Merging image " + str(image_count + 1) + ": " + folder
                
                # Load the current image as a 16 bit numpy array
                current_image = Image.open(os.path.join(datasets_folder,folder,file))
                
                current_image_array = numpy.asarray(current_image).astype('uint16')

                if image_count == 0:
                    # Load the initial image into a numpy 16 bit array
                    base_image_array = current_image_array
                else:
                    
                    # Add the current image to the base image
                    base_image_array = base_image_array + current_image_array
                                    
                image_count = image_count + 1
                
                # Clean up objects in memory
                del current_image
                del current_image_array

    # Get the max pixel value
    pixel_max = base_image_array.max()
      
    # Calculate the divisor for the base image to scale back to 8 bit
    divisor = pixel_max / float(255)
        
    # Calculate back to 8 bit range
    print "Scaling composite image back to 8 bit range (0-255).  Current max pixel value is " + str(pixel_max) + " so will divide by " + str(divisor) + "."
    base_image_array /= divisor
        
    # Convert the array to 8 bit
    print "Converting to 8 bit image"
    base_image_array_8bit = base_image_array.astype('uint8')
        
    # Save the image
    print "Saving image file"
    result_image = Image.fromarray(base_image_array_8bit)
    result_image.save("map.png")
    del result_image
    
    print "Composite image complete"            

def compute_polygon_burn_value(area):
    
    global bounding_box_area, area_range, x_max, x_min, y_max, y_min, resolution
    
    # Compute the area of the bounding box
    bounding_box_area = ((x_max - x_min) * (y_max - y_min))
            
    # ----------------------------------
    # Settings for burn value decay rate
    # ----------------------------------
    
    # Maximum burn value for a polygon
    max_burn_value = 180
    
    # Starting polygon size decay rate begins to be applied
    polygon_area_start_decay = 4 * 2590000  # Roughly four square miles (not really since we are in web mercator but a precise value doesn't matter here) 
        
    # Maximum polygon size relative to the bounding bounding area to have a burn weight
    polygon_area_end_decay = float(bounding_box_area) * float(0.1) # Polygons over 40% of the area of the bounding box won't have any burn value
    
    # Decay rate of burn value (exponential) as polygons increase in size from 0
    decay_rate = 4
        
    burn_value = max_burn_value
    if (area > polygon_area_start_decay and area <= polygon_area_end_decay):
                
        # Range of potential area values
        area_range = polygon_area_end_decay - polygon_area_start_decay
        
        # Calculate percentage of burn value to use
        # We want to weight of the smaller polygons to be much higher so we use
        # the sq. root of the inverse percentage of the area to rapidly unweight 
        #
        area_delta_percentage = float(1) - (float(area) / float(area_range))**(float(1) / float(decay_rate))
        
        # Get the burn value as the inverse percentage of the total
        burn_value = int(area_delta_percentage * float(max_burn_value))
        
        #print str(area) + " of " + str(bounding_box_area) + ".  Percentage: " + str(area_delta_percentage * 100) + "% Burn value:" + str(burn_value)
    elif area > polygon_area_end_decay:
        burn_value = 0
        
    return burn_value

def colorize_composite_image():
    
    # Open the composite image
    image = Image.open("map.png")
    
    # Add an alpha channel
    image = image.convert("RGB")
    
    # Open the color scheme
    scheme = Image.open(os.path.join("scheme", "opencolorado.png"))
    
    print "Loading color scheme.."
    
    # Convert the scheme into a list of clor value tuples
    colors = []
    for row in range(256):
        scheme_pixel = scheme.getpixel((0,255-row))
           
        colors.append(scheme_pixel)
            
    # Load the image pixels for rapid updating
    pixels = image.load()
        
    # Iterate through all of pixels in the composite image
    w,h = image.size
 
    print "Colorizing density map.."
    
    for x in range(w):
        for y in range(h):
            
            # Get the pixel at the location
            pixel = pixels[x,y]
                        
            # Get the color to set the pixel based on the value of the first band
            color = list(colors[pixel[0]])
                                      
            pixels[x,y] = (color[0], color[1], color[2])
 
    print "Saving colorized density map.."
               
    image.save("map_colorized.png")

def generate_world_file():
    
    global resolution, x_min, y_max
        
    world_file = open("map_colorized.pgw", 'wt')
    world_file.write("%0.8f\n" % resolution)  # width of an image pixel in world coordinates
    world_file.write("%0.8f\n" % 0) # Rotation (n/a)
    world_file.write("%0.8f\n" % 0) # Rotation (n/a)
    world_file.write("%0.8f\n" % -resolution) # height of an image pixel in world coordinates (negative)
    world_file.write("%0.8f\n" % x_min) #X coordinate in world coordinates of the center of the top-left pixel
    world_file.write("%0.8f\n" % y_max) #Y coordinate in world coordinates of the center of the top-left pixel
    world_file.close()
    
#Execute main function    
if __name__ == '__main__':
    main()