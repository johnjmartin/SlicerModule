import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy, math

#
# JohnMartin
#

class JohnMartin(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "JohnMartin" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# JohnMartinWidget
#

class JohnMartinWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # threshold value
    #
    self.imageThresholdSliderWidget = ctk.ctkSliderWidget()
    self.imageThresholdSliderWidget.singleStep = 0.1
    self.imageThresholdSliderWidget.minimum = -100
    self.imageThresholdSliderWidget.maximum = 100
    self.imageThresholdSliderWidget.value = 0.5
    self.imageThresholdSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    parametersFormLayout.addRow("Image threshold", self.imageThresholdSliderWidget)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    self.enableScreenshotsFlagCheckBox.checked = 0
    self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    #em selector
    self.emSelector = slicer.qMRMLNodeComboBox()
    self.emSelector.nodeTypes = ['vtkMRMLLinearTransformNode']
    self.emSelector.setMRMLScene(slicer.mrmlScene)
    parametersFormLayout.addRow("EM tool tip transform: ", self.emSelector)

    #opticalSelector

    self.opticalSelector = slicer.qMRMLNodeComboBox()
    self.opticalSelector.nodeTypes = ['vtkMRMLLinearTransformNode']
    self.opticalSelector.setMRMLScene(slicer.mrmlScene)
    parametersFormLayout.addRow("Optical tool tip transform: ", self.opticalSelector)


    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    # self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    # self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.emSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.opticalSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.emSelector.currentNode() and self.opticalSelector.currentNode()

  def onApplyButton(self):
    emTipTransform = self.emSelector.currentNode()
    if emTipTransform == None:
      return
    opTipTransform = self.opticalSelector.currentNode()
    if opTipTransform == None:
      return

    emTipTransform.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, self.onTransformModified)
    opTipTransform.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, self.onTransformModified)

  def onTransformModified(self, caller, event):
    emTipTransform = self.emSelector.currentNode()
    if emTipTransform == None:
      return
    opTipTransform = self.opticalSelector.currentNode()
    if opTipTransform == None:
      return

    emTip_EmTip = [0, 0, 0, 1]
    opTip_OpTip = [0, 0, 0, 1]

    emTipToRasMatrix = vtk.vtkMatrix4x4()
    emTipTransform.GetMatrixTransformToWorld(emTipToRasMatrix)
    emTip_Ras = numpy.array( emTipToRasMatrix.MultiplyFloatPoint(emTip_EmTip) )

    opTipToRasMatrix = vtk.vtkMatrix4x4()
    opTipTransform.GetMatrixTransformToWorld(opTipToRasMatrix)
    opTip_Ras = numpy.array( opTipToRasMatrix.MultiplyFloatPoint(opTip_OpTip) )

    distance = numpy.linalg.norm(emTip_Ras - opTip_Ras)
    print distance

#
# JohnMartinLogic
#


class JohnMartinLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def averageTransformedDistance(self, alphaPoints, betaPoints, alphaToBetaMatrix):
      average = 0
      num = 0

      numberOfPoints = alphaPoints.GetNumberOfPoints()
      bNum = betaPoints.GetNumberOfPoints()

      if numberOfPoints != bNum:
          logging.error('number of points in two lists do not match')
          return -1

      for i in range(numberOfPoints):
          num = num + 1
          a = alphaPoints.GetPoint(i)
          pointA_Alpha = numpy.array(a)
          pointA_Alpha = numpy.append(pointA_Alpha, 1)
          pointA_Beta = alphaToBetaMatrix.MultiplyFloatPoint(pointA_Alpha)
          b = betaPoints.GetPoint(i)
          pointB_Beta = numpy.array(b)
          pointB_Beta = numpy.append(pointB_Beta, 1)
          distance = numpy.linalg.norm(pointA_Beta - pointB_Beta)
          average = average+ (distance-average) / num

      return average

  def rigidRegistration(self, alphaPoints, betaPoints, alphaToBetaMatrix):
    landmarkTransform = vtk.vtkLandmarkTransform()
    landmarkTransform.SetSourceLandmarks(alphaPoints)
    landmarkTransform.SetTargetLandmarks(betaPoints)
    landmarkTransform.SetModeToRigidBody()
    landmarkTransform.Update()
    landmarkTransform.GetMatrix(alphaToBetaMatrix)



class JohnMartinTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def generatePoints(self, numPoints, Scale, Sigma):
    rasFids = slicer.util.getNode('RasPoints')
    if rasFids == None:
        rasFids = slicer.vtkMRMLMarkupsFiducialNode()
        rasFids.SetName('RasPoints')
        slicer.mrmlScene.AddNode(rasFids)
    rasFids.RemoveAllMarkups()

    refFids = slicer.util.getNode('ReferencePoints')
    if refFids == None:
        refFids = slicer.vtkMRMLMarkupsFiducialNode()
        refFids.SetName('ReferencePoints')
        slicer.mrmlScene.AddNode(refFids)
    refFids.RemoveAllMarkups()
    refFids.GetDisplayNode().SetSelectedColor(1,1,0)

    fromNormCoordinates = numpy.random.rand(numPoints, 3)
    noise = numpy.random.normal(0.0, Sigma, numPoints*3)

    for i in range(numPoints):
        x = (fromNormCoordinates[i, 0] - 0.5) * Scale
        y = (fromNormCoordinates[i, 1] - 0.5) * Scale
        z = (fromNormCoordinates[i, 2] - 0.5) * Scale
        rasFids.AddFiducial(x, y, z)
        xx = x+noise[i*3]
        yy = y+noise[i*3+1]
        zz = z+noise[i*3+2]
        refFids.AddFiducial(xx, yy, zz)

  def fiducialsToPoints(self, fiducials, points):
      n = fiducials.GetNumberOfFiducials()
      for i in range(n):
          p = [0,0,0]
          fiducials.GetNthFiducialPosition(i, p)
          points.InsertNextPoint(p[0], p[1], p[2])

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_JohnMartin1()

  def test_JohnMartin1(self):
    # Switch to a layout (24) that contains a Chart View to initiate the construction of the widget and Chart View Node
    lns = slicer.mrmlScene.GetNodesByClass('vtkMRMLLayoutNode')
    lns.InitTraversal()
    ln = lns.GetNextItemAsObject()
    ln.SetViewArrangement(24)

    # Get the Chart View Node
    cvns = slicer.mrmlScene.GetNodesByClass('vtkMRMLChartViewNode')
    cvns.InitTraversal()
    cvn = cvns.GetNextItemAsObject()

    self.delayDisplay("Starting the test")

    referenceToRas = slicer.vtkMRMLLinearTransformNode()
    referenceToRas.SetName('ReferenceToRas')
    slicer.mrmlScene.AddNode(referenceToRas)

    createModelsLogic = slicer.modules.createmodels.logic()
    rasCoordinateModel = createModelsLogic.CreateCoordinate(25, 2)
    rasCoordinateModel.SetName('RasCoordinateModel')
    referenceCoordinateModel = createModelsLogic.CreateCoordinate(20, 2)
    referenceCoordinateModel.SetName('ReferenceCoordinateModel')
    referenceCoordinateModel.SetAndObserveTransformNodeID(referenceToRas.GetID())

    rasCoordinateModel.GetDisplayNode().SetColor(1, 0, 0)
    referenceCoordinateModel.GetDisplayNode().SetColor(0, 0, 1)

    refPoints = vtk.vtkPoints()
    rasPoints = vtk.vtkPoints()

    logic = JohnMartinLogic()

    tre_list = []

    for i in range(10):
        N = 10 + i * 5
        Sigma = 3
        Scale = 100
        self.generatePoints(N, Scale, Sigma)
        rasFids = slicer.util.getNode("RasPoints")
        refFids = slicer.util.getNode('ReferencePoints')
        self.fiducialsToPoints(rasFids, rasPoints)
        self.fiducialsToPoints(refFids, refPoints)
        referenceToRasMatrix = vtk.vtkMatrix4x4()
        logic.rigidRegistration(refPoints, rasPoints, referenceToRasMatrix)
        det = referenceToRasMatrix.Determinant()
        if det < 1e-8:
            print 'Unstable registration. Check input for collinear points.'
            continue
        referenceToRas.SetMatrixTransformToParent(referenceToRasMatrix)
        avgDistance = logic.averageTransformedDistance(refPoints, rasPoints, referenceToRasMatrix)
        print "Avg Distance: " + str(avgDistance)
        targetPoint_Ras = numpy.array([0,0,0,1])
        targetPoint_Reference = referenceToRasMatrix.MultiplyFloatPoint(targetPoint_Ras)
        targetPoint_Reference = numpy.array(targetPoint_Reference)
        tre = numpy.linalg.norm(targetPoint_Ras - targetPoint_Reference)
        tre_list.append(tre)
        print "TRE: " + str(tre)
        print ""

    # Create an Array Node and add some data
    dn = slicer.mrmlScene.AddNode(slicer.vtkMRMLDoubleArrayNode())
    a = dn.GetArray()
    a.SetNumberOfTuples(10)
    x = range(0, 10)
    for i in range(len(x)):
        a.SetComponent(i, 0, (10 + i * 5))
        a.SetComponent(i, 1, tre_list[i])
        a.SetComponent(i, 2, 0)

    cn = slicer.mrmlScene.AddNode(slicer.vtkMRMLChartNode())
    # Add the Array Nodes to the Chart. The first argument is a string used for the legend and to refer to the Array when setting properties.
    cn.AddArray('TRE', dn.GetID())

    # Set a few properties on the Chart. The first argument is a string identifying which Array to assign the property.
    # 'default' is used to assign a property to the Chart itself (as opposed to an Array Node).
    cn.SetProperty('default', 'title', 'Total registration error function')
    cn.SetProperty('default', 'xAxisLabel', 'Points in registration')
    cn.SetProperty('default', 'yAxisLabel', 'TRE')

    # Tell the Chart View which Chart to display
    cvn.SetChartNodeID(cn.GetID())
