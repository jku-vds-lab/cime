import { combineReducers } from 'redux';
import currentTool from "../Ducks/CurrentToolDuck";
import activeStory from "../Ducks/ActiveStoryDuck";
import projectionOpen from "../Ducks/ProjectionOpenDuck";
import highlightedSequence from "../Ducks/HighlightedSequenceDuck";
import dataset from "../Ducks/DatasetDuck";
import openTab from "../Ducks/OpenTabDuck";
import webGLView from "../Ducks/WebGLViewDuck";
import clusterMode from "../Ducks/ClusterModeDuck";
import advancedColoringSelection from "../Ducks/AdvancedColoringSelectionDuck";
import projectionColumns from "../Ducks/ProjectionColumnsDuck";
import displayMode from "../Ducks/DisplayModeDuck";
import lineBrightness from "../Ducks/LineBrightnessDuck";
import activeLine from "../Ducks/ActiveLineDuck";
import stories from "../Ducks/StoriesDuck";
import storyMode from "../Ducks/StoryModeDuck";
import currentAggregation from "../Ducks/AggregationDuck";
import selectedClusters from "../Ducks/SelectedClustersDuck";
import { viewTransform } from "../Ducks/ViewTransformDuck";
import currentClusters from "../Ducks/CurrentClustersDuck";
import projectionParams from "../Ducks/ProjectionParamsDuck";
import checkedShapes from "../Ducks/CheckedShapesDuck";
import projectionWorker from "../Ducks/ProjectionWorkerDuck";
import vectorByShape from "../Ducks/VectorByShapeDuck";
import clusterEdges from "../Ducks/ClusterEdgesDuck";
import selectedVectorByShape from "../Ducks/SelectedVectorByShapeDuck";
import pathLengthRange from '../Ducks/PathLengthRange';
import categoryOptions from '../Ducks/CategoryOptionsDuck';
import channelSize from '../Ducks/ChannelSize';
import globalPointSize from '../Ducks/GlobalPointSizeDuck';
import hoverState from '../Ducks/HoverStateDuck';

export const rootReducer = combineReducers({
  currentTool: currentTool,
  currentAggregation: currentAggregation,
  activeStory: activeStory,
  stories: stories,
  currentClusters: currentClusters,
  openTab: openTab,
  clusterEdges: clusterEdges,
  selectedVectorByShape: selectedVectorByShape,
  vectorByShape: vectorByShape,
  checkedShapes: checkedShapes,
  activeLine: activeLine,
  dataset: dataset,
  highlightedSequence: highlightedSequence,
  viewTransform: viewTransform,
  advancedColoringSelection: advancedColoringSelection,
  storyMode: storyMode,
  projectionColumns: projectionColumns,
  projectionOpen: projectionOpen,
  projectionParams: projectionParams,
  projectionWorker: projectionWorker,
  webGLView: webGLView,
  clusterMode: clusterMode,
  selectedClusters: selectedClusters,
  displayMode: displayMode,
  lineBrightness: lineBrightness,
  pathLengthRange: pathLengthRange,
  categoryOptions: categoryOptions,
  channelSize: channelSize,
  globalPointSize: globalPointSize,
  hoverState: hoverState
})