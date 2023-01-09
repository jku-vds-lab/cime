/**
 * Duck file for the RDKit Settings
 */

const SET_CONTOURLINES = 'ducks/rdkitsettings/SET_CONTOURLINES';
const SET_SCALE = 'ducks/rdkitsettings/SET_SCALE';
const SET_SIGMA = 'ducks/rdkitsettings/SET_SIGMA';
const SET_REFRESH = 'ducks/rdkitsettings/SET_REFRESH';
const SET_SHOW_MCS = 'ducks/rdkitsettings/SET_SHOW_MCS';
const SET_WIDTH = 'ducks/rdkitsettings/SET_WIDTH';
const SET_DO_ALIGNMENT = 'ducks/rdkitsettings/SET_DO_ALIGNMENT';
const SET_DO_DOMAIN = 'ducks/rdkitsettings/SET_DO_DOMAIN';
const SET_SHOW_BARCHART = 'ducks/rdkitsettings/SHOW_BARCHART';
const SET_SHOW_ATTRIBUTES = 'ducks/rdkitsettings/SHOW_ATTRIBUTES';

export const setRDKit_contourLines = (input) => ({
  type: SET_CONTOURLINES,
  input: input,
});

export const setRDKit_scale = (input) => ({
  type: SET_SCALE,
  input: input,
});

export const setRDKit_sigma = (input) => ({
  type: SET_SIGMA,
  input: input,
});

export const setRDKit_refresh = (input) => ({
  type: SET_REFRESH,
  input: input,
});

export const setRDKit_showMCS = (input) => ({
  type: SET_SHOW_MCS,
  input: input,
});

export const setRDKit_width = (input) => ({
  type: SET_WIDTH,
  input: input,
});

export const setRDKit_doAlignment = (input) => ({
  type: SET_DO_ALIGNMENT,
  input: input,
});

export const setRDKit_colorDomain = (input) => ({
  type: SET_DO_DOMAIN,
  input: input,
});

export const setRDKit_showBarChart = (input) => ({
  type: SET_SHOW_BARCHART,
  input: input,
});

export const setRDKit_showAttributes = (input) => ({
  type: SET_SHOW_ATTRIBUTES,
  input: input,
});

const initialState: RDKitSettingsType = {
  contourLines: 10,
  scale: -1,
  sigma: 0,
  refresh: 0,
  showMCS: false,
  width: 250,
  doAlignment: true,
  domain: {
    type: 0,
    value: null,
    deadzone: 0,
    thresholds: [],
  },
  showBarChart: false,
  showAttributes: false,
};

export type RDKitSettingsType = {
  contourLines: number;
  scale: number;
  sigma: number;
  refresh: number;
  showMCS: boolean;
  width: number;
  doAlignment: boolean;
  domain: {
    type: number;
    value: number[];
    deadzone: number;
    thresholds: number[];
  };
  showBarChart: boolean;
  showAttributes: boolean;
};

const rdkitSettings = (state = initialState, action): RDKitSettingsType => {
  switch (action.type) {
    case SET_CONTOURLINES:
      return { ...state, contourLines: action.input };
    case SET_SCALE:
      return { ...state, scale: action.input };
    case SET_SIGMA:
      return { ...state, sigma: action.input };
    case SET_REFRESH:
      return { ...state, refresh: action.input };
    case SET_SHOW_MCS:
      return { ...state, showMCS: action.input };
    case SET_WIDTH:
      return { ...state, width: action.input };
    case SET_DO_ALIGNMENT:
      return { ...state, doAlignment: action.input };
    case SET_DO_DOMAIN:
      return {
        ...state,
        domain: { ...state.domain, ...action.input },
      };
    case SET_SHOW_BARCHART:
      return {
        ...state,
        showBarChart: action.input,
      };
    case SET_SHOW_ATTRIBUTES:
      return {
        ...state,
        showAttributes: action.input,
      };
    default:
      return state;
  }
};

export default rdkitSettings;
