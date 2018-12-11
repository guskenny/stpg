The aim of the files here is to make it relatively easy for a people to pull out 
a project on one of the OR windows servers (eg ifap-dev-or, quijibo, babylon-bp) 
and compile a project in Visual Studio that links agains local copies of relevant
libraries such as boost, coin, cplex, gurobi.


How to set up a NEW PROJECT:
* Set up a project in Visual Studio as normal
* Save & close the visual studio project and edit the .vcxproj file (XML file) in  
  a text editor.
* add the following XML after   <ImportGroup Label="ExtensionSettings"></ImportGroup>
  <ImportGroup Label="PropertySheets">
  	<Import Project="vc_config\$(COMPUTERNAME)_$(PlatformToolset).props" Condition="exists('vc_config\$(COMPUTERNAME)_$(PlatformToolset).props')" Label="LocalAppDataPlatform" />
  </ImportGroup>
  This should suck in the correct vc_config\IFAP-DEV-OR_v120.props property sheet or similar
* Add vc_config directory as an "external" directory to the directory containing the .vcxproj 
  file. This can be done using TortoiseSVN -> Properties -> New... -> Externals
  Then "New..." with URL: https://svnserv.csiro.au/svn/mos/shared/vc_config/
* Use the user macros defined in the property sheets such as $(boostLib) to refer to the libraries
  to be linked into the project. (All of the include directories should automatically be added)
  
  
How to configure a NEW COMPUTER:
Most things are being done via Macros. To see what macros are defined in Visual Studio,
right click on a project, bring up properties, select <Edit...> for any propoperty and 
click the "Macros>>>" button. This should show for example COMPUTERNAME

* Copy an existing property sheet (eg IFAP-DEV-OR_v120.props) and rename it as COMPUTERNAME_v120.props
* The copied sheet should automatically be loaded into a project (provided it is configured
  as above)
* Check the "Property Manager" tab (next to "Solution Explorer" in Visual Studio) where it  
  should appear under "project name" --> Release | x64 (or any of the other configurations).
* The values can be edited by right clicking on the sheet and selecting "Properties"
* Most of the options are controlled by the "User Macros" part
* When in doubt check the .props file in a text/XML editor


