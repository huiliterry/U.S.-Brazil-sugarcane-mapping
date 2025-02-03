/******S2 image collection******/
var S2_cutCldSlw = function(start_date,end_date,boundary,CLD_PRB_THRESH,NIR_DRK_THRESH,CLD_PRJ_DIST,BUFFER){
  // print('CLD_PRB_THRESH',CLD_PRB_THRESH);
  // Function to add NDVI, time, and constant variables to Sentinel imagery.
  var addVariables = function(image) {
    // Compute time in fractional years since the epoch.
    var date = ee.String(image.get('system:index'));
    // var days = ee.Date(date.slice(0,4).cat('-').cat(date.slice(4,6)).cat('-').cat(date.slice(6,8))).format('DDD');
    var year = date.slice(0,4);
    var month = date.slice(4,6);
    var dateOfMonth = date.slice(6,8);
    //generating day number in the imgcollection
    var days = ee.Number.parse(ee.Date(year.cat('-').cat(month).cat('-').cat(dateOfMonth)).format('DDD'))
                        .add((ee.Number.parse(year).subtract(ee.Number.parse(startYear))).multiply(365))
                        .subtract(initialDays);

    // Return the image with the added bands.
    return image
    // Add an NDVI band.
    .addBands(image.normalizedDifference(['B8', 'B4']).float().rename('NDVI'))
    // Add an GCVI band.
    .addBands(image.select('B8').divide(image.select('B3')).subtract(ee.Image(1)).float().rename('GCVI'))
    // Add an NDMI band.
    .addBands(image.normalizedDifference(['B8', 'B12']).float().rename('NDMI'))
    // Add an MSI band.
    .addBands(image.select('B11').divide(image.select('B8')).float().rename('MSI'))
    // Add an NDWI band.
    .addBands(image.normalizedDifference(['B3','B8']).rename(['NDWI']))
    // edit band names.
    .select(band)
    .set('system:day_start',ee.Number.parse(days));
  };
  /************************************************************/
  // Remove clouds, add variables and filter to the area of interest.
  // Cloud components
  var add_cloud_bands = function(img){
      //Get s2cloudless image, subset the probability band.
      var cld_prb = ee.Image(img.get('s2cloudless')).select('probability');
  
      //Condition s2cloudless by the probability threshold value.
      var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds');
  
      //Add the cloud probability layer and cloud mask as image bands.
      return img.addBands(ee.Image([cld_prb, is_cloud]));
  };
  //Cloud shadow components
  var add_shadow_bands = function(img){
      //Identify water pixels from the SCL band.
      //var not_water = img.select('SCL').neq(6)
  
      //Identify dark NIR pixels that are not water (potential cloud shadow pixels).
      var SR_BAND_SCALE = 1e4;
      var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).rename('dark_pixels');//.multiply(not_water)
  
      //Determine the direction to project cloud shadow from clouds (assumes UTM projection).
      var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));
  
      //Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
      var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
          .reproject({'crs': img.select(0).projection(), 'scale': 100})
          .select('distance')
          .mask()
          .rename('cloud_transform'));
  
      //Identify the intersection of dark pixels with cloud shadow projection.
      var shadows = cld_proj.multiply(dark_pixels).rename('shadows');
      //Add dark pixels, cloud projection, and identified shadows as image bands.
      return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]));
  };
  //Final cloud-shadow mask
  var add_cld_shdw_mask = function(img){
      //Add cloud component bands.
      var img_cloud = add_cloud_bands(img);
      //Add cloud shadow component bands.
      var img_cloud_shadow = add_shadow_bands(img_cloud);
      //Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
      var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0);
      //Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
      //20 m scale is for speed, and assumes clouds don't require 10 m precision.
      is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
          .reproject({'crs': img.select([0]).projection(), 'scale': 20})
          .rename('cloudmask'));
      //Add the final cloud-shadow mask to the image.
      return img_cloud_shadow.addBands(is_cld_shdw);
  };
  
  var apply_cld_shdw_mask = function(img){
      //Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
      var not_cld_shdw = img.select('cloudmask').not();
  
      //Subset reflectance bands and update their masks, return the result.
      return img.select('B.*').updateMask(not_cld_shdw);
  };
  //mask the water pixel
  var apply_scl_water_mask = function(img){
    var scl = img.select('SCL');
    var wantedPixels = scl.neq(6);
    var targetPixels = scl.eq(4).or(scl.eq(5));
    return img.updateMask(wantedPixels);
  };
  // Import and filter S2 SR.
  var s2_sr_col = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')//S2_SR_HARMONIZED
          .filterDate(start_date, end_date)
          .filterBounds(boundary)
          // .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER))
          .filter(ee.Filter.eq('GENERAL_QUALITY','PASSED'));
  // Import and filter s2cloudless.
  var s2_cloudless_col = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
          .filterBounds(boundary)
          .filterDate(start_date, end_date)  ;  
  var s2_sr_cld_col_eval = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
          'primary': s2_sr_col,
          'secondary': s2_cloudless_col,
          'condition': ee.Filter.equals({
              'leftField': 'system:index',
              'rightField': 'system:index'
          })
      }));
  //finding the initial days of the year    
  var initialDate = ee.String(s2_sr_col.first().get('system:index'));
  var startYear = initialDate.slice(0,4);
  var month = initialDate.slice(4,6);
  var dateOfMonth = initialDate.slice(6,8);
  var initialDays = ee.Number.parse(ee.Date(startYear.cat('-').cat(month).cat('-').cat(dateOfMonth)).format('DDD'));
  // print('initialDays',initialDays)
  //creating final imgcollection
  var s2_no_cld_shdw =  s2_sr_cld_col_eval//s2_sr_cld_col_eval //s2_sr_col
                        .map(add_cld_shdw_mask)
                        .map(apply_cld_shdw_mask)
                        // .map(apply_scl_water_mask)
                        .map(function(img){
                          return img.clip(boundary);//.divide(10000);
                        })
                        .map(addVariables);
  // print('s2_no_cld_shdw',s2_no_cld_shdw);
  return s2_no_cld_shdw;
};
var compositeS2  = function(parameters){
  parameters = parameters || {};
  var startDate = parameters.startDate;
  var endDate = parameters.endDate;
  // var S2_CLOUD_FILTER = parameters.S2_CLOUD_FILTER;
  var S2L89_CLD_PRB_THRESH = parameters.S2L89_CLD_PRB_THRESH;
  var S2_NIR_DRK_THRESH = parameters.S2_NIR_DRK_THRESH;
  var S2_CLD_PRJ_DIST = parameters.S2_CLD_PRJ_DIST;
  var S2_BUFFER = parameters.S2_BUFFER;
  // var L89_CLOUD_FILTER = parameters.L89_CLOUD_FILTER;
  var area_boundary = parameters.area_boundary;
  var intervalValue = parameters.intervalValue;
  
  var preProsS2 = S2_cutCldSlw(startDate,endDate,area_boundary,S2L89_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER);
  
  var merge15Interval = intervalYear(startDate,endDate,area_boundary,preProsS2,intervalValue);
  return merge15Interval;
};
/******cosine regression******/
//calculate harmonic coefficients
var ResCoesHarm = function(order,imgCollection){
  var dependentSeries = ee.List(band);
  // The number of cycles per year to model.
  // Make a list of harmonic frequencies to model.  
  // These also serve as band name suffixes.
  var harmonicFrequencies = ee.List.sequence(1, order);
  // Function to get a sequence of band names for harmonic terms.
  var getNames = function(base, list) {
    return ee.List(list).map(function(i) { 
      return ee.String(base).cat(ee.Number(i).int());
    });
  };
  // Construct lists of names for the harmonic terms.
  var cosNames = getNames('cos_', harmonicFrequencies);
  // Independent variables.
  var independents = ee.List(['constant']).cat(cosNames);
  // print('independents',independents);

  var addConstant = function(image) {
    return image.addBands(ee.Image(1));
  };
  var basicFrequency = 0.5;
  // Function to add a time band.
  var addTime = function(image) {
    // Compute time in fractional years since the epoch.
    var days = ee.String(image.get('system:day_start'));
    var timeRadians = ee.Image(
      ee.Number.parse(days).divide(365).multiply(basicFrequency*2*Math.PI)
    );
    return image.addBands(timeRadians.rename('t').float());
  };
  // addTime(s2_no_cld_shdw.first());
  
  var addHarmonics = function(freqs) {
    return function(image) {
      // Make an image of frequencies.
      var frequencies = ee.Image.constant(freqs);
      // This band should represent time in radians.
      var time = ee.Image(image).select('t');
      // Get the sin terms.
      var cosines = time.multiply(frequencies).cos().rename(cosNames);
      return image.addBands(cosines);
    };
  };
  //print('s2_no_cld_shdw',s2_no_cld_shdw)
  var harmonicS2 = imgCollection
    .map(addConstant)
    .map(addTime)
    .map(addHarmonics(harmonicFrequencies));
  // print('harmonicS2',harmonicS2);
  //amplitudes
  var amplitudes = ee.ImageCollection(dependentSeries.map(function(dependent){
    return harmonicS2.select(independents.add(dependent))
    .reduce(ee.Reducer.robustLinearRegression(independents.length(), 1))
    .select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents])
    .slice(1)
    .set('spectral',dependent)
    .rename([
              ee.String(dependent).cat(ee.String("_cos_1")),ee.String(dependent).cat(ee.String("_cos_2")),
              ee.String(dependent).cat(ee.String("_cos_3")),ee.String(dependent).cat(ee.String("_cos_4")),
              ee.String(dependent).cat(ee.String("_cos_5")),ee.String(dependent).cat(ee.String("_cos_6")),
              ee.String(dependent).cat(ee.String("_cos_7")),ee.String(dependent).cat(ee.String("_cos_8")),
              ee.String(dependent).cat(ee.String("_cos_9")),ee.String(dependent).cat(ee.String("_cos_10")),
              ee.String(dependent).cat(ee.String("_cos_11")),ee.String(dependent).cat(ee.String("_cos_12"))
             ]);
  }));
  return amplitudes;
}; 
var coefficientsImg = function(startDate,endDate,testBounds,S2L89_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order){
  var preProsS2 = compositeS2L89({startDate:startDate,
                                  endDate:endDate,
                                  S2L89_CLD_PRB_THRESH:S2L89_CLD_PRB_THRESH,
                                  S2_NIR_DRK_THRESH:S2_NIR_DRK_THRESH,
                                  S2_CLD_PRJ_DIST:S2_CLD_PRJ_DIST,
                                  S2_BUFFER:S2_BUFFER,
                                  area_boundary:testBounds,
                                  intervalValue:intervalValue
                  });
  //generating harmonic coefficients
  var coefficients = ResCoesHarm(order,preProsS2);
  return coefficients;
};
/******classifier training******/
var bandsName;
var trainingSampleConstr = function(startDateTrain,endDateTrain,startDateTrain_NS,endDateTrain_NS,FloridaLabels,LouisianaLabels,order,S2L89_CLD_PRB_THRESH){//,forestThreshold,Maxmin
   var S2_NIR_DRK_THRESH = 0.25;
  var S2_CLD_PRJ_DIST = 1;
  var S2_BUFFER = 20; 
  //preparing the coefficients, training dataset, classifer
  //samples geometry
  var sugarcane_Florida = ee.FeatureCollection('TIGER/2018/Counties').filter(ee.Filter.eq('COUNTYNS','00295761')).geometry();
  var sugarcane_Louisiana = ee.FeatureCollection('TIGER/2018/Counties').filter(ee.Filter.eq('COUNTYNS','00558065')).geometry();
  var nonSugarcaneRegions = ee.FeatureCollection('TIGER/2018/Counties')
  .filter(ee.Filter.or(
    ee.Filter.eq('COUNTYNS','00758521'),
    ee.Filter.eq('COUNTYNS','00758553'),
    ee.Filter.eq('COUNTYNS','00758556')
  )).union().geometry();
  
  //samples coefficients
  var FloridacoefficientsImgs = coefficientsImg(startDateTrain,endDateTrain,sugarcane_Florida,S2L89_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order).toBands();
  print('FloridacoefficientsImgs',FloridacoefficientsImgs);
  var LouisianacoefficientsImgs = coefficientsImg(startDateTrain,endDateTrain,sugarcane_Louisiana,S2L89_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order).toBands();
  print('LouisianacoefficientsImgs',LouisianacoefficientsImgs);
  
  bandsName = FloridacoefficientsImgs.bandNames();
  // generating the training dataset using the above sample.
  var Florida_trainingSample = FloridacoefficientsImgs.sampleRegions({
    collection: FloridaLabels,//.merge(randomPoints_nonS),
    properties: ['label'],
    scale: 10,
    tileScale:4
  });
  
  var Louisiana_trainingSample = LouisianacoefficientsImgs.sampleRegions({
    collection: LouisianaLabels,//.merge(randomPoints_nonS),
    properties: ['label'],
    scale: 10,
    tileScale:4
  });  
  var sugarcaneSamples = Louisiana_trainingSample.merge(Florida_trainingSample);
  
  var randomPoints_nonS = ee.FeatureCollection.randomPoints({
        region: nonSugarcaneRegions, points: sugarcaneSamples.size(), seed: 0, maxError: 1//.divide(1.5).round().multiply(2) 
    })
    .map(function(feature){return feature.set('label',0);});

  
  var coefficientsImgs_nonS = coefficientsImg(startDateTrain_NS,endDateTrain_NS,nonSugarcaneRegions,S2L89_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order).toBands();
 
  var trainingSample_nonS = coefficientsImgs_nonS.sampleRegions({
    collection: randomPoints_nonS,
    properties: ['label'],
    scale: 10,
    tileScale:4
  });
  return sugarcaneSamples.merge(trainingSample_nonS);
};
var classifierConstr = function(trainingSamples,bandsName){
  var classifier = ee.Classifier.smileRandomForest(100).train({
                      features: trainingSamples,
                      classProperty: 'label',
                      inputProperties: bandsName
                    });
  return classifier;
};
