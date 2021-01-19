// show/hide plots and table
$(document).ready(function(){
$('#dtBin').click(function(){
  $('#DtBin').toggle();
});

$('#intConfig1').click(function(){
    $('#DtBin').slideUp();
  $('#IntConfig1').toggle();
});

$('#intConfig2').click(function(){
  $('#IntConfig2').toggle();
});

$('#vwIntraMap').click(function(){
  $('#IntraMap').toggle();
	//$('#dSwitch').toggleClass('hide');
  //$('#Map1').toggleClass('col-sm-6 col-sm-4');
  //$('#Map2').toggleClass('col-sm-6 col-sm-8');
});

$('#vwMap1').click(function(){
  $('#IntraMap').slideUp();
  $('#Map1').toggle();
});

$('#vwScaf').click(function(){
  $('#IntraMap').slideUp();
  $('#VwScaf').toggle();
});


$('#vwJointMap').click(function(){
  $('#VwJointMap').toggle();
});

$('#vwIntData').click(function(){
  $('#VwIntData').toggle();
    $('#VwJointMap').slideUp();
});

$('#vwJointMap2').click(function(){
  $('#VwJointMap').slideUp();
   $('#VwIntData').slideUp();
  $('#VwJointMap2').toggle();
});

$('#scafConfig').click(function(){
  $('#IntConfig2').slideUp();
  $('#ScafConfig').toggle();
});

  // reset the app when the button pressed
  $('#reset').click(function(){
  history.go(0);
});


// tooltips
$("#dir1").attr('title', 'Scaffolding direction');
$("#mapSrc").attr('title', 'Data source');
$("#exitZoom").attr('title', 'Exit zoom');
$("#clearBrush").attr('title', 'Clear cut line');
$("#cut2").attr('title', 'Cut sequence');
});

