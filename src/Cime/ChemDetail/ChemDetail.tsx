import * as React from "react";
import "./chem.scss";
import {
  Autocomplete,
  Box,
  Button,
  Checkbox,
  createFilterOptions,
  FormControl,
  FormControlLabel,
  FormGroup,
  Grid,
  IconButton,
  Input,
  InputLabel,
  Paper,
  Popover,
  Slider,
  Switch,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import { trackPromise } from "react-promise-tracker";
import { connect, ConnectedProps } from "react-redux";
import RefreshIcon from "@mui/icons-material/Refresh";
import SettingsIcon from "@mui/icons-material/Settings";
import InfoIcon from "@mui/icons-material/Info";
import FilterListIcon from "@mui/icons-material/FilterList";
import { isFunction } from "lodash";
import { useCancellablePromise } from "projection-space-explorer";
import AddCircleOutlineIcon from "@mui/icons-material/AddCircleOutline";
import DeleteIcon from "@mui/icons-material/Delete";
import {
  selectVectors,
  setHoverState,
  WindowMode,
} from "projection-space-explorer";
import { CIMEBackendFromEnv } from "../../Backend/CIMEBackend";
import { LoadingIndicatorView } from "../../Overrides/DatasetTabPanel";
import { AppState } from "../../State/Store";
import {
  setRDKit_contourLines,
  setRDKit_doAlignment,
  setRDKit_refresh,
  setRDKit_scale,
  setRDKit_showMCS,
  setRDKit_sigma,
  setRDKit_width,
} from "../../State/RDKitSettingsDuck";

/**
 * Chem Legend, implemented
 */

const mapStateToProps_Chem = (state: AppState) => ({
  dataset: state.dataset,
  hoverSettings: state.hoverSettings,
  rdkitSettings: state.rdkitSettings,
  columns: state.dataset?.columns,
});
const mapDispatchToProps_Chem = (dispatch: any) => ({
  setCurrentAggregation: (samples: number[]) =>
    dispatch(selectVectors(samples, false)),
  setHoverstate: (state: any, updater: any) =>
    dispatch(setHoverState(state, updater)),
});
const connector_Chem = connect(mapStateToProps_Chem, mapDispatchToProps_Chem);

/**
 * Type that holds the props we declared above in mapStateToProps and mapDispatchToProps
 */
type PropsFromRedux_Chem = ConnectedProps<typeof connector_Chem>;

type Props_Chem_Parent = PropsFromRedux_Chem & {
  selection: any[];
  aggregate: boolean;
  mcs_only?: boolean;
  diff?: boolean;
  selection_ref?: any[];
};

export const ChemLegendParent = connector_Chem(function (
  props: Props_Chem_Parent
) {
  const { cancellablePromise, cancelPromises } = useCancellablePromise();
  if (props.mcs_only) {
    const [mcsComp, setMcsComp] = React.useState(<div>loading...</div>);

    let smiles_col = get_smiles_col(props.columns);

    React.useEffect(() => {
      cancelPromises();

      if (smiles_col in props.columns) {
        const controller = new AbortController();
        let my_fetch: Promise<any>;

        if (props.diff && props.selection_ref) {
          const smilesA = props.selection.map((row) => row[smiles_col]);
          const smilesB = props.selection_ref.map((row) => row[smiles_col]);
          my_fetch = CIMEBackendFromEnv.getDifferenceHighlight(
            smilesA,
            smilesB,
            controller
          );
        } else {
          const formData = new FormData();
          props.selection.every((row) => {
            formData.append("ids", row["cime_id"]);
            return true;
          });
          my_fetch = CIMEBackendFromEnv.getMCSFromSmilesList(
            formData,
            controller
          );
        }

        trackPromise(
          cancellablePromise(
            my_fetch.then((x) => {
              if (x.length > 100) {
                // check if it is actually long enogh to be an img
                setMcsComp(() => (
                  <div
                    style={{
                      width: 200,
                      height: 200,
                      backgroundSize: "contain",
                      backgroundPosition: "center",
                      backgroundRepeat: "no-repeat",
                      backgroundImage: `url('data:image/jpg;base64,${x}')`,
                    }}
                  ></div>
                ));
              } else {
                setMcsComp(() => <div>{x}</div>);
              }
            }),
            controller
          )
        );
      }
    }, [props.selection, props.selection_ref, props.mcs_only]);

    if (smiles_col in props.columns) {
      return <div>{mcsComp}</div>;
    }
    return <div>No SMILES column found.</div>;
  }

  const [settingsOpen, setSettingsOpen] = React.useState(false);
  const [repList, setRepList] = React.useState(["Common Substructure"]);
  const [chemComponents, setChemComponents] = React.useState([0]);
  const [chemComponentsCurrentRep, setChemComponentsCurrentRep] =
    React.useState(["Common Substructure"]);

  const loadRepList = function (refresh = false) {
    if (refresh || repList.length <= 1) {
      const loading_area = "global_loading_indicator";
      const controller = new AbortController();
      trackPromise(
        cancellablePromise(
          CIMEBackendFromEnv.getRepresentationList(
            refresh,
            props.dataset.info.path,
            controller
          ).then((x) => {
            if (x["rep_list"].length > 0) {
              let rep_list = [...x["rep_list"]];
              rep_list.splice(0, 0, "Common Substructure");
              setRepList(rep_list);
            }
          }),
          controller
        ),
        loading_area
      );
    }
  };

  const addComp = function () {
    let comps = [...chemComponents];
    comps.push(Math.max(...comps) + 1);
    setChemComponents(comps);

    let compsCR = [...chemComponentsCurrentRep];
    compsCR.push("Common Substructure");
    setChemComponentsCurrentRep(compsCR);
  };

  React.useEffect(() => {
    cancelPromises();
    if (props.aggregate) {
      loadRepList();
    }
  }, []);

  const removeComponent = (id) => {
    let comps = [...chemComponents];
    let compsCR = [...chemComponentsCurrentRep];
    const index = comps.indexOf(id);
    if (index > -1) {
      comps.splice(index, 1);
      compsCR.splice(index, 1);
    }
    setChemComponents(comps);
    setChemComponentsCurrentRep(compsCR);
  };

  const setCurrentRep = (value, id) => {
    if (repList.includes(value)) {
      let compsCR = [...chemComponentsCurrentRep];
      const index = chemComponents.indexOf(id);
      compsCR[index] = value;
      setChemComponentsCurrentRep(compsCR);
    }
  };

  const anchorRef = React.useRef<any>();
  const chemRef = React.useRef<any>();

  if (props.aggregate) {
    return (
      <Box className={"ParentChem"}>
        {props.aggregate && (
          <Box paddingLeft={2} paddingRight={2}>
            <Tooltip title="Summary Settings">
              <Button
                style={{ color: "gray" }}
                ref={anchorRef}
                onClick={() => setSettingsOpen(true)}
              >
                <SettingsIcon></SettingsIcon>&nbsp; Settings
              </Button>
            </Tooltip>
            <SettingsPopover
              open={settingsOpen}
              setOpen={setSettingsOpen}
              anchorEl={anchorRef.current}
              refreshRepList={() => {
                loadRepList(true);
              }}
            ></SettingsPopover>
            <Tooltip title="Add Component">
              <Button style={{ color: "gray" }} onClick={() => addComp()}>
                <AddCircleOutlineIcon></AddCircleOutlineIcon>&nbsp; Add View
              </Button>
            </Tooltip>
          </Box>
        )}
        <div ref={chemRef} className={"chemComponents"}>
          {chemComponents.length > 1 && (
            <div
              style={{
                width: (props.rdkitSettings.width + 20) * chemComponents.length,
              }}
            >
              {chemComponents.map((x, i) => {
                return (
                  <div
                    key={x}
                    style={{
                      width: props.rdkitSettings.width + 20,
                      float: "left",
                    }}
                  >
                    <ChemLegend
                      chemRef={chemRef}
                      setCurrentRep={(value) => setCurrentRep(value, x)}
                      currentRep={chemComponentsCurrentRep[i]}
                      removeComponent={() => removeComponent(x)}
                      id={x}
                      rep_list={repList}
                      selection={props.selection}
                      aggregate={props.aggregate}
                    ></ChemLegend>
                  </div>
                );
              })}
            </div>
          )}
          {chemComponents.length <= 1 && (
            <div>
              <div
                style={{ minWidth: props.rdkitSettings.width }}
                key={chemComponents[0]}
              >
                <ChemLegend
                  chemRef={chemRef}
                  setCurrentRep={(value) =>
                    setCurrentRep(value, chemComponents[0])
                  }
                  currentRep={chemComponentsCurrentRep[0]}
                  id={chemComponents[0]}
                  rep_list={repList}
                  selection={props.selection}
                  aggregate={props.aggregate}
                ></ChemLegend>
              </div>
            </div>
          )}
        </div>
      </Box>
    );
  } else {
    return (
      <ChemLegend
        id={-1}
        rep_list={repList}
        selection={props.selection}
        aggregate={props.aggregate}
      ></ChemLegend>
    );
  }
});

type Props_Chem = PropsFromRedux_Chem & {
  selection: any;
  columns: any;
  aggregate: boolean;
  rep_list: string[];
  id: any;
  removeComponent?: any;
  currentRep?: string;
  setCurrentRep?: any;
  chemRef?;
};

const loading_area = "chemlegend_loading_area";
const UPDATER = "chemdetail";
const ChemLegend = connector_Chem(
  class extends React.Component<Props_Chem, { checkedList: boolean[] }> {
    anchorRef: any;

    constructor(props) {
      super(props);
      this.state = {
        checkedList: [],
      };
    }

    render() {
      const handleMouseEnter = (i) => {
        let hover_item = null;
        if (i >= 0) {
          hover_item = this.props.selection[i];
        }
        this.props.setHoverstate(hover_item, UPDATER);
      };

      const handleMouseOut = () => {
        let hover_item = null;
        this.props.setHoverstate(hover_item, UPDATER);
      };

      const setCheckedList = (value) => {
        const set_val = isFunction(value)
          ? value(this.state.checkedList)
          : value;
        this.setState({ ...this.state, checkedList: set_val });
      };

      const handle_filter = () => {
        const filter_instances = this.props.selection.filter(
          (x, i) => this.state.checkedList[i]
        );
        if (filter_instances.length > 0) {
          setCheckedList([]);
          this.props.setCurrentAggregation(
            filter_instances.map((e) => e.__meta__.meshIndex)
          );
        } else {
          alert(
            "Please, select at least one Compound in the Summary View to filter."
          );
        }
      };

      if (this.props.aggregate) {
        return (
          <div className={"ParentImg"}>
            <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
              <RepresentationList
                value={this.props.currentRep}
                onChange={this.props.setCurrentRep}
                rep_list={this.props.rep_list}
                hoverSettings={this.props.hoverSettings}
              />
            </Box>

            <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
              <Button
                size="small"
                variant="outlined"
                onClick={() => {
                  handle_filter();
                }}
              >
                <FilterListIcon fontSize={"small"} />
                &nbsp;Confirm Selection
              </Button>
              {this.props.removeComponent && (
                <IconButton onClick={this.props.removeComponent}>
                  <DeleteIcon></DeleteIcon>
                </IconButton>
              )}
            </Box>
            <LoadingIndicatorView area={loading_area + this.props.id} />
            <ImageView
              chemRef={this.props.chemRef}
              id={this.props.id}
              setCheckedList={setCheckedList}
              selection={this.props.selection}
              columns={this.props.columns}
              aggregate={this.props.aggregate}
              current_rep={this.props.currentRep}
              handleMouseEnter={handleMouseEnter}
              handleMouseOut={handleMouseOut}
            />
          </div>
        );
      }

      return (
        <div>
          <ImageView
            id={this.props.id}
            selection={this.props.selection}
            columns={this.props.columns}
            aggregate={this.props.aggregate}
          />
        </div>
      );
    }
  }
);

function get_smiles_col(columns) {
  let smiles_col = "SMILES";
  // TODO: find by meta_data -> how to handle multiple smiles columns? for now: just take first column that contains "smiles"
  if (!(smiles_col in columns)) {
    let col_names = Object.keys(columns);
    for (const key in col_names) {
      let col = col_names[key];
      if (col.toLowerCase().includes("smiles")) {
        smiles_col = col;
        break;
      }
    }
  }
  return smiles_col;
}

const mapStateToProps_Img = (state: AppState) => ({
  datasetState: state.dataset,
  hoverState: state.hoverState,
  rdkitSettings: state.rdkitSettings,
});
const mapDispatchToProps_Img = (dispatch) => ({});
const connector_Img = connect(mapStateToProps_Img, mapDispatchToProps_Img);

/**
 * Type that holds the props we declared above in mapStateToProps and mapDispatchToProps
 */
type PropsFromRedux_Img = ConnectedProps<typeof connector_Img>;

type Props_Img = PropsFromRedux_Img & {
  selection;
  columns;
  aggregate;
  handleMouseEnter?;
  handleMouseOut?;
  current_rep?;
  setCheckedList?;
  id;
  chemRef?;
};

function addHighlight(element) {
  if (element && element.style) {
    element.style["border"] = "solid black 1px";
  }
}

function removeHighlight(element) {
  if (element && element.style) {
    element.style["border"] = "solid white 1px";
  }
}

const ImageView = connector_Img(function ({
  chemRef,
  id,
  hoverState,
  selection,
  columns,
  aggregate,
  handleMouseEnter,
  handleMouseOut,
  current_rep,
  setCheckedList,
  rdkitSettings,
  datasetState,
}: Props_Img) {
  const [comp, setComp] = React.useState(<div></div>);

  const ref = React.useRef<any>();
  const { cancellablePromise, cancelPromises } = useCancellablePromise();

  React.useEffect(() => {
    cancelPromises(); // cancel all unresolved promises
  }, [selection, current_rep]);

  React.useEffect(() => {
    if (setCheckedList) setCheckedList([]);

    let smiles_col = get_smiles_col(columns);

    const onUpdateItem = (i, val) => {
      setCheckedList((checkedList) => {
        const list = checkedList.map((item, j) => {
          if (j === i) {
            return val;
          } else {
            return item;
          }
        });
        return list;
      });
    };

    // TODO: find by meta_data -> how to handle multiple smiles columns?
    // for (const col_name in columns) {
    //     let col = columns[col_name];
    //     if(col.metaInformation.imgSmiles)
    //         smiles_col = col_name;
    // }
    if (smiles_col in columns) {
      setComp(<div></div>);
      if (selection.length > 0) {
        if (aggregate) {
          const formData = new FormData();
          formData.append("current_rep", current_rep);
          selection.forEach((row) => {
            formData.append("ids", row["cime_id"]);
          });
          formData.append(
            "contourLines",
            rdkitSettings.contourLines?.toString()
          );
          formData.append("scale", rdkitSettings.scale?.toString());
          formData.append("sigma", rdkitSettings.sigma?.toString());
          formData.append("showMCS", rdkitSettings.showMCS?.toString());
          formData.append("width", rdkitSettings.width?.toString());
          formData.append("doAlignment", rdkitSettings.doAlignment?.toString());
          formData.append("filename", datasetState.info.path);

          const controller = new AbortController();
          trackPromise(
            cancellablePromise(
              CIMEBackendFromEnv.getStructuresFromSmilesList(
                formData,
                controller
              ),
              controller
            ).then((x) => {
              const img_lst = x["img_lst"].map((base64, i) => {
                setCheckedList((checkedList) => {
                  let cpy_checked_list = [...checkedList];
                  if (cpy_checked_list.length <= i) {
                    cpy_checked_list.push(false);
                  }
                  return cpy_checked_list;
                });
                return (
                  <Grid className={"legend_multiple"} key={i} item>
                    <FormControlLabel
                      labelPlacement="bottom"
                      control={
                        <Checkbox
                          color="primary"
                          onChange={(event) => {
                            onUpdateItem(i, event.target.checked);
                          }}
                        />
                      }
                      label={
                        <img
                          src={"data:image/jpeg;base64," + base64}
                          alt=""
                          onMouseEnter={() => {
                            handleMouseEnter(i);
                          }}
                          onMouseOver={() => {
                            handleMouseEnter(i);
                          }}
                          onMouseLeave={() => {
                            handleMouseOut();
                          }}
                        />
                      }
                    />
                    <Typography style={{ paddingLeft: 5 }} variant="subtitle2">
                      ID: {selection[i]["cime_id"]}
                    </Typography>
                  </Grid>
                );
              }); //key={props.selection[i][smiles_col]} --> gives error because sometimes smiles ocure twice
              //<div dangerouslySetInnerHTML={{ __html: img_lst.join("") }} />
              setComp(img_lst);
            }),
            loading_area + id
          );
        } else {
          let row = selection[0];
          const controller = new AbortController();
          cancellablePromise(
            CIMEBackendFromEnv.getStructureFromSmiles(
              datasetState.info.path,
              row[smiles_col],
              false,
              controller
            ),
            controller
          )
            .then((x) => {
              setComp(
                <div>
                  <img
                    className={"legend_single"}
                    src={"data:image/jpeg;base64," + x}
                    alt=""
                  />
                  <Typography style={{ paddingLeft: 5 }} variant="subtitle2">
                    ID: {row["cime_id"]}
                  </Typography>
                </div>
              );
            })
            .catch((error) => console.log(error));
        }
      } else {
        setComp(<div>No Selection</div>);
      }
    } else {
      setComp(<div>No SMILES column found</div>);
    }
  }, [selection]);

  React.useEffect(() => {
    if (aggregate) {
      let smiles_col = get_smiles_col(columns);

      if (smiles_col in columns) {
        let imgList = ref.current.childNodes;
        if (selection.length === imgList.length) {
          ref.current.style.display = "none";

          const formData = new FormData();
          formData.append("current_rep", current_rep);
          selection.forEach((row) => {
            formData.append("ids", row["cime_id"]);
          });
          formData.append(
            "contourLines",
            rdkitSettings.contourLines?.toString()
          );
          formData.append("scale", rdkitSettings.scale?.toString());
          formData.append("sigma", rdkitSettings.sigma?.toString());
          formData.append("showMCS", rdkitSettings.showMCS?.toString());
          formData.append("width", rdkitSettings.width?.toString());
          formData.append("doAlignment", rdkitSettings.doAlignment?.toString());
          formData.append("filename", datasetState.info.path);

          const controller = new AbortController();
          trackPromise(
            cancellablePromise(
              CIMEBackendFromEnv.getStructuresFromSmilesList(
                formData,
                controller
              ),
              controller
            ).then((x) => {
              x["img_lst"].forEach((base64, i) => {
                const cur_img = imgList[i].getElementsByTagName("img")[0];
                cur_img.src = "data:image/jpeg;base64," + base64;
                // cur_img.width = props.rdkitSettings.width;
                // cur_img.height = props.rdkitSettings.width;
              });
              ref.current.style.display = "flex";
            }),
            loading_area + id
          );
        }
      }
    }
  }, [current_rep, rdkitSettings.refresh]);

  React.useEffect(() => {
    if (aggregate) {
      //@ts-ignore
      let container = chemRef?.current;
      let imgContainer = container.getElementsByClassName("chem-grid")[0];
      //@ts-ignore
      let imgList = imgContainer.childNodes;
      if (hoverState && hoverState.data) {
        const idx = selection.findIndex(
          (x) =>
            x &&
            x["__meta__"] &&
            hoverState.data["__meta__"] &&
            x["__meta__"]["meshIndex"] ==
            hoverState.data["__meta__"]["meshIndex"]
        );
        if (idx >= 0 && imgList.length > 0) {
          for (const i in imgList) {
            const img_div = imgList[i];
            removeHighlight(img_div);
          }
          addHighlight(imgList[idx]);

          if (hoverState.updater != UPDATER) {
            if (container && imgList[idx]) {
              //@ts-ignore
              container.scrollTop =
                imgList[idx].offsetTop - container.offsetTop;
              // imgList[idx]?.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'start' }) // this seems to be buggy sometimes
            }
          }
        }
      } else {
        for (const i in imgList) {
          const img_div = imgList[i];
          removeHighlight(img_div);
        }
      }
    }
  }, [hoverState.data, hoverState.updater]);

  return (
    <div className={"chemContainer"}>
      <Grid ref={ref} className={"chem-grid"} container>
        {comp}
      </Grid>
    </div>
  );
});

interface ValLabelProps {
  children: React.ReactElement;
  open: boolean;
  value: number;
}
function ValueLabelComponent(props: ValLabelProps) {
  const { children, open, value } = props;

  return (
    <Tooltip open={open} enterTouchDelay={0} placement="top" title={value}>
      {children}
    </Tooltip>
  );
}

const mapStateToProps_settings = (state: AppState) => ({
  rdkitSettings: state.rdkitSettings,
});
const mapDispatchToProps_settings = (dispatch) => ({
  setContourLines: (input) => dispatch(setRDKit_contourLines(input)),
  setScale: (input) => dispatch(setRDKit_scale(input)),
  setSigma: (input) => dispatch(setRDKit_sigma(input)),
  setShowMCS: (input) => dispatch(setRDKit_showMCS(input)),
  setWidth: (input) => dispatch(setRDKit_width(input)),
  setRefresh: (input) => dispatch(setRDKit_refresh(input)),
  setDoAlignment: (input) => dispatch(setRDKit_doAlignment(input)),
});
const connector_settings = connect(
  mapStateToProps_settings,
  mapDispatchToProps_settings
);

type PropsFromRedux_Settings = ConnectedProps<typeof connector_settings>;

type SettingsPopoverProps = PropsFromRedux_Settings & {
  open: boolean;
  setOpen: any;
  anchorEl: any;
  refreshRepList: any;
};

const SettingsPopover = connector_settings(function ({
  open,
  setOpen,
  anchorEl,
  refreshRepList,
  rdkitSettings,
  setContourLines,
  setScale,
  setSigma,
  setShowMCS,
  setWidth,
  setRefresh,
  setDoAlignment,
}: SettingsPopoverProps) {
  return (
    <Popover
      disablePortal={true}
      id={"dialog to open"}
      open={open}
      anchorEl={anchorEl}
      onClose={() => setOpen(() => false)}
      anchorOrigin={{
        vertical: "bottom",
        horizontal: "right",
      }}
      transformOrigin={{
        vertical: "bottom",
        horizontal: "left",
      }}
    >
      <div>
        <Paper style={{ padding: 10, minWidth: 300 }}>
          <Grid container spacing={2} direction='column' justifyContent='stretch'
            sx={{
              '& .MuiTextField-root': { width: '200px' },
            }}>
            <Grid item>
              <Button
                size="small"
                variant="outlined"
                aria-label={"Refresh Representation List"}
                onClick={() => refreshRepList(true)}
              >
                <RefreshIcon />
                Refresh Representation List
              </Button>
            </Grid>

            <Grid item>
              <Typography variant="subtitle2" gutterBottom>
                RDKit Settings
              </Typography>
            </Grid>

            {/* https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/Draw/SimilarityMaps.py */}
            {/* how many contour lines should be drawn [0;inf] */}

            <Grid item container alignItems='center'>
              <TextField
                label="Contour Lines"
                type="number"
                InputLabelProps={{
                  shrink: true
                }}
                value={rdkitSettings.contourLines}
                onChange={(event) => {
                  let val = parseInt(event.target.value);
                  if (isNaN(val)) setContourLines(event.target.value);
                  else setContourLines(Math.max(val, 0));
                }}

              />

              <div style={{ flexGrow: 1 }}></div>
              <Tooltip title="Number of Contour Lines [0; &infin;] &isin; &#8469;">
                <InfoIcon fontSize="small"></InfoIcon>
              </Tooltip>
            </Grid>

            {/* scale tells the programm about how to scale the weights [-1;inf]; default is -1 which means that it is inherted by the algorithm*/}
            <Grid item container alignItems='center'>
              <TextField
                label="Scale"
                type="number"
                InputLabelProps={{
                  shrink: true
                }}
                value={rdkitSettings.scale}
                onChange={(event) => {
                  let val = parseFloat(event.target.value);
                  if (isNaN(val)) setScale(event.target.value);
                  else setScale(Math.max(val, -1));
                }}

              />

              <div style={{ flexGrow: 1 }}></div>
              <Tooltip title="Weight Scale [-1; &infin;] &isin; &#8477;"><InfoIcon fontSize="small"></InfoIcon></Tooltip>
            </Grid>

            {/* sigma is for gaussian ]0;~0.2?]; default 0 means that the algorithm infers the value from the weights */}
            <Grid item container alignItems='center'>
              <TextField
                label="Sigma"
                type="number"
                InputLabelProps={{
                  shrink: true
                }}
                inputProps={{
                  step: 0.1
                }}
                value={rdkitSettings.sigma}
                onChange={(event) => {
                  let val = parseFloat(event.target.value);
                  if (isNaN(val)) setSigma(event.target.value);
                  else setSigma(Math.max(val, 0));
                }}

              />

              <div style={{ flexGrow: 1 }}></div>
              <Tooltip title="Sigma for Gaussian ]0; &infin;] &isin; &#8477;. Default of 0 signals the algorithm to infer the value.">
                <InfoIcon fontSize="small"></InfoIcon>
              </Tooltip>
            </Grid>


            <Grid item>
              <FormControlLabel
                control={
                  <Switch
                    color="primary"
                    checked={rdkitSettings.showMCS}
                    onChange={(_, value) => {
                      setShowMCS(value);
                    }}
                  />
                }
                label="Show MCS"
              />
            </Grid>

            <Grid item>
              <FormControlLabel
                control={
                  <Switch
                    color="primary"
                    checked={rdkitSettings.doAlignment}
                    onChange={(_, value) => {
                      setDoAlignment(value);
                    }}
                  />
                }
                label="Align Structure"
              />
            </Grid>

            <Grid item>
              <Typography style={{ paddingTop: 10 }} gutterBottom>
                Image Width
              </Typography>

              <Slider
                value={rdkitSettings.width}
                onChange={(event, new_val) => {
                  setWidth(new_val);
                }}
                min={50}
                max={500}
                step={10}
              />
            </Grid>

            <Grid item>
              <Button
                style={{ marginTop: 3, maxWidth: 150 }}
                size="small"
                variant="outlined"
                onClick={() => {
                  setRefresh((rdkitSettings.refresh += 1));
                }}
              >
                Apply Settings
              </Button>
            </Grid>
          </Grid>
        </Paper>
      </div>
    </Popover>
  );
});

const RepresentationList = (props) => {
  const options = props.rep_list.map((rep) => {
    let split = rep.split("_");
    const inputVal = split.pop();
    let group = split.join("_");
    group = group.replace("atom.dprop.", "");
    group = group.replace("atom.dprop", "");
    return {
      group: group,
      value: rep,
      inputValue: inputVal,
    };
  });

  const filterOptions = createFilterOptions({
    stringify: (option: any) => {
      return option.value;
    },
  });


  return (
    <Autocomplete
      size={"small"}
      className={props.className}
      filterOptions={filterOptions}
      onChange={(event, newValue) => {
        if (newValue) props.onChange(newValue.value);
      }}
      disablePortal={props.hoverSettings.windowMode == WindowMode.Extern}
      options={options.sort((a, b) => -b.group.localeCompare(a.group))}
      groupBy={(option: any) => option.group}
      getOptionLabel={(option: any) => option.inputValue}
      style={{ maxWidth: 300 }}
      defaultValue={options[0]}
      renderInput={(params) => (
        <TextField
          {...params}
          label="Choose Representation"
          variant="outlined"
        />
      )}
      isOptionEqualToValue={(option, value) => option.value == value.value}
    />
  );
};
//getOptionSelected={(option:any, value) => {return option.value == value.value;}}
