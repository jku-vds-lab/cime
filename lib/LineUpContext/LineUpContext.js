var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
import { connect } from "react-redux";
import * as LineUpJS from "lineupjs";
import "./LineUpContext.scss";
import { createSelectionDesc, Column, ERenderMode, renderMissingDOM, StringColumn, } from "lineupjs";
import { CIMEBackendFromEnv } from "../Backend/CIMEBackend";
import * as _ from "lodash";
import React from "react";
import { ACluster, DiscreteMapping, EXCLUDED_COLUMNS, MyWindowPortal, PrebuiltFeatures, selectVectors, setDetailVisibility, setHoverState, ShallowSet, AStorytelling, } from "projection-space-explorer";
import { TestColumn } from "./LineUpClasses/TestColumn";
import { setLineUpInput_lineup } from "../State/LineUpInputDuck";
import * as d3v5 from "d3v5";
import isEqual from "lodash/isEqual";
/**
 * Declares a function which maps application state to component properties (by name)
 *
 * @param state The whole state of the application (contains a field for each duck!)
 */
var mapStateToProps = function (state) {
    var _a, _b;
    var activeMultiple = state.multiples.multiples.entities[state.multiples.active];
    return {
        dataset: state.dataset,
        lineUpInput: state.lineUpInput,
        lineUpInput_data: (_a = state.dataset) === null || _a === void 0 ? void 0 : _a.vectors,
        lineUpInput_columns: (_b = state.dataset) === null || _b === void 0 ? void 0 : _b.columns,
        currentAggregation: state.currentAggregation,
        activeStory: AStorytelling.getActive(state.stories),
        pointColorScale: activeMultiple === null || activeMultiple === void 0 ? void 0 : activeMultiple.attributes.pointColorScale,
        channelColor: activeMultiple === null || activeMultiple === void 0 ? void 0 : activeMultiple.attributes.channelColor,
        detailView: state.detailView,
        // splitRef: state.splitRef
        //hoverState: state.hoverState
    };
};
/**
 * Declares a function which maps dispatch events to component properties (by name)
 *
 * @param dispatch The generic dispatch function declared in redux
 */
var mapDispatchToProps = function (dispatch) { return ({
    setCurrentAggregation: function (samples) {
        return dispatch(selectVectors(samples));
    },
    setLineUpInput_visibility: function (visibility) {
        return dispatch(setDetailVisibility(visibility));
    },
    setLineUpInput_lineup: function (input) { return dispatch(setLineUpInput_lineup(input)); },
    setHoverstate: function (state, updater) { return dispatch(setHoverState(state, updater)); },
}); };
/**
 * Factory method which is declared here so we can get a static type in 'ConnectedProps'
 */
var connector = connect(mapStateToProps, mapDispatchToProps);
function arrayEquals(a, b) {
    return (Array.isArray(a) &&
        Array.isArray(b) &&
        a.length === b.length &&
        a.every(function (val, index) { return val === b[index]; }));
}
// let lineup = null;
var UPDATER = "lineup";
var UNIQUE_ID = "unique_ID";
/**
 * Our component definition, by declaring our props with 'Props' we have static types for each of our property
 */
export var LineUpContext = connector(function (_a) {
    var _b;
    var dataset = _a.dataset, lineUpInput = _a.lineUpInput, lineUpInput_data = _a.lineUpInput_data, lineUpInput_columns = _a.lineUpInput_columns, currentAggregation = _a.currentAggregation, channelColor = _a.channelColor, setCurrentAggregation = _a.setCurrentAggregation, setLineUpInput_lineup = _a.setLineUpInput_lineup, setLineUpInput_visibility = _a.setLineUpInput_visibility, onFilter = _a.onFilter, activeStory = _a.activeStory, pointColorScale = _a.pointColorScale, setHoverstate = _a.setHoverstate, detailView = _a.detailView;
    // In case we have no input, dont render at all
    if (!lineUpInput || !lineUpInput_data || !detailView.open) {
        //splitRef?.current?.setSizes([100, 0])
        return null;
    }
    var lineup_ref = React.useRef();
    var debouncedHighlight = React.useCallback(_.debounce(function (hover_item) { return setHoverstate(hover_item, UPDATER); }, 200), []);
    var preprocess_lineup_data = function (data) {
        if (activeStory)
            ACluster.deriveVectorLabelsFromClusters(data, Object.values(activeStory.clusters.entities));
        var lineup_data = new Array();
        var columns = {};
        data.forEach(function (element) {
            // if(element[PrebuiltFeatures.ClusterLabel].length <= 0){
            //     element[PrebuiltFeatures.ClusterLabel] = [-1];
            // }
            var row = {};
            for (var i in lineUpInput_columns) {
                var col = lineUpInput_columns[i];
                if (!EXCLUDED_COLUMNS.includes(i) &&
                    (Object.keys(col.metaInformation).length <= 0 ||
                        !col.metaInformation.noLineUp)) {
                    if (Object.keys(col.metaInformation).length > 0 &&
                        col.metaInformation.lineUpGroup) {
                        // if(col.metaInformation.lineUpGroup.endsWith(""))
                        var split = col.metaInformation.lineUpGroup.split("_");
                        if (split.length <= 1) {
                            // if the string is separated with and underscore, only the first part of the string is considered as the group. the second part of the string determines a sub value of this group
                            if (Object.keys(row).includes(col.metaInformation.lineUpGroup)) {
                                row[col.metaInformation.lineUpGroup].push(element[i]);
                            }
                            else {
                                row[col.metaInformation.lineUpGroup] = [element[i]];
                                columns[col.metaInformation.lineUpGroup] = col;
                            }
                        }
                        else {
                            var group_name = split[0];
                            var var_name = split[1];
                            if (Object.keys(row).includes(group_name)) {
                                if (Object.keys(row[group_name]).includes(var_name)) {
                                    row[group_name][var_name].push(element[i]);
                                }
                                else {
                                    row[group_name][var_name] = [element[i]];
                                }
                            }
                            else {
                                row[group_name] = {};
                                row[group_name][var_name] = [element[i]];
                            }
                            // update column metaInformation
                            if (Object.keys(columns).includes(group_name)) {
                                columns[group_name].metaInformation.globalMin = Math.min(columns[group_name].metaInformation.globalMin, element[i]);
                                columns[group_name].metaInformation.globalMax = Math.max(columns[group_name].metaInformation.globalMax, element[i]);
                            }
                            else {
                                columns[group_name] = col;
                                columns[group_name].metaInformation.customLineChart = true;
                                columns[group_name].metaInformation.globalMin = element[i];
                                columns[group_name].metaInformation.globalMax = element[i];
                            }
                        }
                    }
                    else {
                        row[i] = element[i];
                        columns[i] = col;
                    }
                }
            }
            row[PrebuiltFeatures.ClusterLabel] =
                element[PrebuiltFeatures.ClusterLabel].toString();
            row[UNIQUE_ID] = element["__meta__"]["meshIndex"];
            lineup_data.push(row);
            // console.log(element)
            // let row = Object.assign({}, element)
            // row[PrebuiltFeatures.ClusterLabel] = element[PrebuiltFeatures.ClusterLabel].toString();
            // row[UNIQUE_ID] = element["__meta__"]["view"]["meshIndex"];
            // lineup_data.push(row);
        });
        return [lineup_data, columns];
    };
    var clear_automatic_filters = function (lineUpInput, filter) {
        if (filter) {
            var _loop_1 = function (key) {
                var lineup = lineUpInput.lineup;
                var ranking = lineup.data.getFirstRanking();
                if (key === "selection") {
                    var filter_col = ranking.children.find(function (x) {
                        return x.desc.column == UNIQUE_ID;
                    });
                    filter_col === null || filter_col === void 0 ? void 0 : filter_col.clearFilter();
                }
                else {
                    var filter_col = ranking.children.find(function (x) {
                        return x.desc.column == key;
                    });
                    filter_col === null || filter_col === void 0 ? void 0 : filter_col.clearFilter();
                }
            };
            for (var key in filter) {
                _loop_1(key);
            }
        }
    };
    var get_lineup_dump = function (lineUpInput) {
        if (lineUpInput.lineup) {
            clear_automatic_filters(lineUpInput, lineUpInput.filter);
            var dump = lineUpInput.lineup.dump();
            return dump;
        }
        return null;
    };
    React.useEffect(function () {
        // if(lineUpInput.dump){
        //     try {
        //         const json_parsed = JSON.parse(lineUpInput.dump)
        //         const restored = fromDumpFile(json_parsed)
        //         console.log(restored);
        //         const builder = buildLineup(lineUpInput.columns, restored.dat).restore(restored.dump);
        //         // const builder = LineUpJS.builder(restored.data).restore(restored.dump);
        //         lineup?.destroy();
        //         lineup = builder.build(lineup_ref.current);
        //         return;
        //     } catch (error) {
        //         console.log(error);
        //     }
        // }
        var _a;
        var _b = preprocess_lineup_data(lineUpInput_data), lineup_data = _b[0], columns = _b[1];
        var builder = buildLineup(columns, lineup_data, pointColorScale, channelColor);
        var dump = get_lineup_dump(lineUpInput);
        (_a = lineUpInput.lineup) === null || _a === void 0 ? void 0 : _a.destroy();
        var lineup = builder.buildTaggle(lineup_ref.current);
        if (dump) {
            lineup.restore(dump);
        }
        var ranking = lineup.data.getFirstRanking();
        // add selection checkbox column
        var selection_col = ranking.children.find(function (column) {
            return column.desc.type === "selection";
        });
        if (!selection_col) {
            selection_col = lineup.data.create(createSelectionDesc("Selections"));
            if (selection_col) {
                ranking.insert(selection_col, 1);
            }
        }
        // // make lineup filter interact with the scatter plot view
        // ranking.on('orderChanged.custom', (previous, current, previousGroups, currentGroups, dirtyReason) => {
        //     if (dirtyReason.indexOf('filter') === -1) {
        //         return;
        //     }
        //     const onRankingChanged = (current) => {
        //         for (let i=0; i < lineUpInput.data.length; i++) {
        //             lineUpInput.data[i].view.lineUpFiltered = !current.includes(i);
        //         }
        //         onFilter()
        //     }
        //     onRankingChanged(current)
        // });
        // make lineup selection interact with the scatter plot view
        lineup.on("selectionChanged", function (currentSelection_lineup) {
            // if(currentSelection_lineup.length == 0) return; // selectionChanged is called during creation of lineup, before the current aggregation was set; therefore, it would always set the current aggregation to nothing because in the lineup table nothing was selected yet
            var currentSelection_scatter = lineUpInput_data
                .map(function (x, i) {
                if (x.__meta__.selected)
                    return i;
            })
                .filter(function (x) { return x !== undefined; });
            if (!arrayEquals(currentSelection_lineup, currentSelection_scatter)) {
                // need to check, if the current lineup selection is already the current aggregation
                var agg_1 = new Array();
                currentSelection_lineup.forEach(function (index) {
                    agg_1.push(lineUpInput_data[index].__meta__.meshIndex);
                });
                setCurrentAggregation(agg_1);
            }
        });
        lineup.on("highlightChanged", function (idx) {
            var hover_item;
            if (idx >= 0) {
                hover_item = lineUpInput_data[idx];
            }
            debouncedHighlight(hover_item);
        });
        // update lineup when smiles_column width changes
        if (smiles_structure_columns.length > 0 ||
            custom_chart_columns.length > 0) {
            var custom_chart_cols = ranking.children.filter(function (x) {
                return custom_chart_columns.includes(x.label);
            });
            for (var i in custom_chart_cols) {
                var custom_chart_col = custom_chart_cols[i];
                custom_chart_col.on("widthChanged", function (prev, current) {
                    lineup.update();
                });
            }
            var lineup_smiles_cols = ranking.children.filter(function (x) { return x instanceof StructureImageColumn; });
            var _loop_2 = function (i) {
                var lineup_smiles_col = lineup_smiles_cols[i];
                lineup_smiles_col.on("widthChanged", function (prev, current) {
                    lineup.update();
                });
                // custom filter adapted from michael
                var filterChanged = function (prev, cur) {
                    if ((prev === null || prev === void 0 ? void 0 : prev.filter) !== (cur === null || cur === void 0 ? void 0 : cur.filter)) {
                        // only update, if it is a new filter
                        var filter_1 = typeof (cur === null || cur === void 0 ? void 0 : cur.filter) === "string" ? cur === null || cur === void 0 ? void 0 : cur.filter : null; // only allow string filters -> no regex (TODO: remove regex checkbox)
                        if (lineup_smiles_col && filter_1) {
                            CIMEBackendFromEnv.getSubstructureCount(
                            // @ts-ignore
                            lineUpInput_data.map(function (d) { return d[lineup_smiles_col.desc.column]; }), filter_1)
                                .then(function (matches) {
                                var validSmiles = matches
                                    .filter(function (_a) {
                                    var smiles = _a[0], count = _a[1];
                                    return count > 0;
                                })
                                    .map(function (_a) {
                                    var smiles = _a[0], count = _a[1];
                                    return smiles;
                                });
                                lineup_smiles_col.setFilter({
                                    filter: filter_1,
                                    valid: new Set(validSmiles),
                                    filterMissing: cur.filterMissing,
                                });
                            })
                                .catch(function (e) {
                                lineup_smiles_col.setFilter(null);
                            });
                        }
                    }
                };
                lineup_smiles_col.on(StringColumn.EVENT_FILTER_CHANGED, filterChanged);
            };
            for (var i in lineup_smiles_cols) {
                _loop_2(i);
            }
        }
        setLineUpInput_lineup(lineup);
    }, [
        lineUpInput_data,
        lineUpInput_columns,
        activeStory,
        activeStory === null || activeStory === void 0 ? void 0 : activeStory.clusters,
        (_b = activeStory === null || activeStory === void 0 ? void 0 : activeStory.clusters) === null || _b === void 0 ? void 0 : _b.ids.length,
        lineUpInput.update,
    ]);
    // React.useEffect(() => { //TODO: not working...
    //     // update lineup, if current storybook (current cluster) changed
    //     if(lineUpInput.lineup){
    //         const data_provider = lineUpInput.lineup.data;
    //         let lineup_data = preprocess_lineup_data(lineUpInput_data);
    //         console.log("setdata")
    //         data_provider.setData(lineup_data);
    //         const ranking = lineUpInput.lineup.data.getFirstRanking();
    //         const my_col_builder = LineUpJS.buildCategoricalColumn(PrebuiltFeatures.ClusterLabel);
    //         console.log(lineup_data)
    //          ranking.insert(lineUpInput.lineup.data.create(my_col_builder.build(lineup_data)));
    //         ranking.insert(lineUpInput.lineup.data.create(my_col_builder.build([]])));
    //         // const ranking = lineUpInput.lineup.data.getFirstRanking();
    //         // // let cluster_col = ranking.columns.find(x => x.desc.column == PrebuiltFeatures.ClusterLabel);
    //         // // const my_desc = cluster_col.desc;
    //         // // my_desc.categories = ["test"]
    //         // // const my_col = new CategoricalColumn(cluster_col.id, cluster_col.desc)
    //         // const my_col_builder = LineUpJS.buildCategoricalColumn(PrebuiltFeatures.ClusterLabel);
    //         // // console.log()
    //         // ranking.insert(lineUpInput.lineup.data.create(my_col_builder.build(lineup_data))); //asSet(',')
    //         // // data_provider.setData(lineup_data)
    //         // // lineUpInput.lineup.update();
    //         // // lineUpInput.lineup?.setDataProvider(data_provider);
    //         // lineUpInput.lineup.restore(lineUpInput.lineup.dump())
    //         // // console.log(cluster_col.dump())
    //         // // console.log(lineUpInput.lineup.dump())
    //     }
    // }, [activeStory, activeStory?.clusters, activeStory?.clusters?.length]);
    // this effect is allways executed after the component is rendered when currentAggregation changed
    React.useEffect(function () {
        if (lineUpInput.lineup != null) {
            // select those instances that are also selected in the scatter plot view
            if (currentAggregation.aggregation &&
                currentAggregation.aggregation.length > 0) {
                var currentSelection_scatter = lineUpInput_data
                    .map(function (x, i) {
                    if (x.__meta__.selected)
                        return i;
                })
                    .filter(function (x) { return x !== undefined; });
                lineUpInput.lineup.setSelection(currentSelection_scatter);
                //console.log("set selection");
                //console.log(lineUpInput);
                //console.log(currentSelection_scatter);
                // const lineup_idx = lineup.renderer?.rankings[0]?.findNearest(currentSelection_scatter);
                // lineup.renderer?.rankings[0]?.scrollIntoView(lineup_idx);
                // set the grouping to selection checkboxes -> uncomment if this should be automatically if something changes
                // const ranking = lineup.data.getFirstRanking();
                // let selection_col = ranking.children.find(x => x.label == "Selection Checkboxes");
                // ranking.groupBy(selection_col, -1) // remove grouping first
                // ranking.groupBy(selection_col);
            }
            else {
                lineUpInput.lineup.setSelection([]);
            }
        }
    }, [lineUpInput.lineup, currentAggregation]);
    React.useEffect(function () {
        if (lineUpInput.lineup && lineUpInput.lineup.data) {
            var ranking = lineUpInput.lineup.data.getFirstRanking();
            clear_automatic_filters(lineUpInput, lineUpInput.previousfilter);
            if (lineUpInput.filter) {
                var _loop_3 = function (key) {
                    var cur_filter = lineUpInput.filter[key];
                    if (key === "reset" && cur_filter) {
                        ranking.clearFilters();
                    }
                    else if (key === "selection") {
                        var filter_col = ranking.children.find(function (x) {
                            return x.desc.column == UNIQUE_ID;
                        });
                        var regex_str_1 = "";
                        lineUpInput.filter[key].forEach(function (element) {
                            regex_str_1 += "|";
                            regex_str_1 += dataset.vectors[element].__meta__.meshIndex;
                        });
                        regex_str_1 = regex_str_1.substr(1); // remove the leading "|"
                        var my_regex = new RegExp("^(".concat(regex_str_1, ")$"), "i"); // i modifier says that it's not case sensitive; ^ means start of string; $ means end of string
                        filter_col === null || filter_col === void 0 ? void 0 : filter_col.setFilter({
                            filter: my_regex,
                            filterMissing: true,
                        });
                    }
                    else {
                        var filter_col = ranking.children.find(function (x) {
                            return x.desc.column == key;
                        });
                        var my_regex = new RegExp("^(.+,)?".concat(cur_filter, "(,.+)?$"), "i"); // i modifier says that it's not case sensitive; ^ means start of string; $ means end of string
                        filter_col === null || filter_col === void 0 ? void 0 : filter_col.setFilter({
                            filter: my_regex,
                            filterMissing: true,
                        });
                    }
                };
                for (var key in lineUpInput.filter) {
                    _loop_3(key);
                }
            }
        }
    }, [lineUpInput.lineup, lineUpInput.filter]);
    //https://github.com/lineupjs/lineup_app/blob/master/src/export.ts
    return false ? (React.createElement(MyWindowPortal, { onClose: function () {
            var _a;
            (_a = lineUpInput.lineup) === null || _a === void 0 ? void 0 : _a.destroy();
            setLineUpInput_visibility(false);
        } },
        React.createElement("div", { ref: lineup_ref, id: "lineup_view" }))) : (React.createElement("div", { className: "LineUpParent" },
        React.createElement("div", { style: {
                clear: "both",
                position: "absolute",
                top: "1px",
                bottom: 0,
                left: 0,
                right: 0,
                padding: 0,
            }, ref: lineup_ref, id: "lineup_view" })));
});
var WIDTH_HEIGHT_RATIO = 2;
var smiles_structure_columns = new Array();
var custom_chart_columns = new Array();
function myDynamicHeight(data, ranking) {
    return { defaultHeight: 25, height: function () { return 25; }, padding: function () { return 0; } };
    if (smiles_structure_columns.length > 0) {
        var cols = ranking.children.filter(function (x) {
            return smiles_structure_columns.includes(x.label) ||
                custom_chart_columns.includes(x.label);
        });
        if (!cols || cols.length == 0)
            return { defaultHeight: 25, height: function () { return 25; }, padding: function () { return 0; } };
        var col_heights = cols.map(function (x) {
            if (custom_chart_columns.includes(x.label))
                return x.getWidth() / WIDTH_HEIGHT_RATIO; // for chart, the width should be bigger than the height
            return x.getWidth(); // for images it is square
        });
        var col_height_1 = Math.max(Math.max.apply(Math, col_heights), 25); //col.getWidth();
        var height = function (item) {
            return col_height_1;
        };
        var padding = function (item) {
            return 0;
        };
        return { defaultHeight: col_height_1, height: height, padding: padding };
    }
    return { defaultHeight: 25, height: function () { return 25; }, padding: function () { return 0; } };
}
// const base_color = undefined;
var base_color = "#c1c1c1";
// const base_color = "#1f77b4";
function buildLineup(cols, data, pointColorScale, channelColor) {
    // console.log(channelColor) //TODO: update lineup colorscale, if sth changes; TODO: do this for all columns, not just groupLabel
    var groupLabel_cat_color;
    if (channelColor &&
        pointColorScale &&
        channelColor.key === PrebuiltFeatures.ClusterLabel) {
        // console.log(pointColorScale)
        var groupLabel_mapping_1 = new DiscreteMapping(pointColorScale, new ShallowSet(data.map(function (vector) { return vector[PrebuiltFeatures.ClusterLabel]; })));
        //console.log(groupLabel_mapping)
        groupLabel_cat_color = groupLabel_mapping_1.values
            .filter(function (cat) { return cat && cat !== ""; })
            .map(function (cat) {
            return { name: cat, color: groupLabel_mapping_1.map(cat).hex };
        });
    }
    var builder = LineUpJS.builder(data);
    for (var i in cols) {
        var col = cols[i];
        var show = true; //!(typeof col.metaInformation.hideLineUp !== 'undefined' && col.metaInformation.hideLineUp); // hide column if "hideLineUp" is specified -> there is a lineup bug with that option
        // if (!EXCLUDED_COLUMNS.includes(i) && (Object.keys(col.metaInformation).length <= 0 || !col.metaInformation.noLineUp)) { // only if there is a "noLineUp" modifier at this column or thix column is excluded, we don't do anything
        if (col.metaInformation.imgSmiles) {
            var smiles_col = "Structure: " + i;
            smiles_structure_columns.push(smiles_col);
            builder.column(LineUpJS.buildColumn("mySmilesStructureColumn", i)
                .label(smiles_col)
                .renderer("mySmilesStructureRenderer", "mySmilesStructureRenderer")
                .width(50)
                .build([]));
            builder.column(LineUpJS.buildStringColumn(i)
                .width(50)
                .custom("visible", show)
                .color(base_color));
        }
        else if (col.metaInformation.customLineChart) {
            // builder.column(LineUpJS.buildNumberColumn(i).label(i).asMap().renderer("myLineChartRenderer", "myLineChartRenderer").width(50).build([]));
            builder.column(LineUpJS.buildColumn("myLineChartColumn", i)
                .label(i)
                .custom("min", col.metaInformation.globalMin)
                .custom("max", col.metaInformation.globalMax)
                .renderer("myLineChartRenderer", "myLineChartRenderer")
                .width(150)
                .build([]));
            custom_chart_columns.push(i);
        }
        else if (i == PrebuiltFeatures.ClusterLabel) {
            var clust_col = LineUpJS.buildCategoricalColumn(i, groupLabel_cat_color)
                .custom("visible", show)
                .width(70); // .asSet(',')
            builder.column(clust_col);
        }
        else {
            builder.deriveColumns(i);
        }
    }
    // builder.deriveColumns([]);
    builder.column(LineUpJS.buildStringColumn("Annotations").editable().color(base_color));
    builder.column(LineUpJS.buildStringColumn(UNIQUE_ID).width(50).color(base_color)); // we need this to be able to filter by all indices; this ID corresponds to the mesh index
    builder.defaultRanking(true);
    // builder.deriveColors();
    builder.registerRenderer("mySmilesStructureRenderer", new MySmilesStructureRenderer());
    builder.registerRenderer("myLineChartRenderer", new MyLineChartRenderer());
    // builder.registerRenderer("myBarCellRenderer", new BarCellRenderer(true));
    builder.registerColumnType("mySmilesStructureColumn", StructureImageColumn);
    builder.registerColumnType("myLineChartColumn", TestColumn);
    builder.sidePanel(true, true); // collapse side panel by default
    builder.livePreviews({
        filter: false,
    });
    builder.dynamicHeight(myDynamicHeight);
    builder.animated(false);
    return builder;
}
var StructureImageColumn = /** @class */ (function (_super) {
    __extends(StructureImageColumn, _super);
    function StructureImageColumn() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.structureFilter = null;
        return _this;
    }
    StructureImageColumn.prototype.filter = function (row) {
        if (!this.isFiltered()) {
            return true;
        }
        return this.structureFilter.valid.has(this.getLabel(row));
    };
    StructureImageColumn.prototype.isFiltered = function () {
        var _a;
        return this.structureFilter != null && ((_a = this.structureFilter.valid) === null || _a === void 0 ? void 0 : _a.size) > 0;
    };
    StructureImageColumn.prototype.getFilter = function () {
        return this.structureFilter;
    };
    StructureImageColumn.prototype.setFilter = function (filter) {
        if (isEqual(filter, this.structureFilter)) {
            return;
        }
        this.fire([
            StringColumn.EVENT_FILTER_CHANGED,
            Column.EVENT_DIRTY_VALUES,
            Column.EVENT_DIRTY,
        ], this.structureFilter, (this.structureFilter = filter));
    };
    return StructureImageColumn;
}(StringColumn));
export { StructureImageColumn };
var MyLineChartRenderer = /** @class */ (function () {
    function MyLineChartRenderer() {
        this.title = "Line Chart";
    }
    MyLineChartRenderer.prototype.canRender = function (col, mode) {
        // return col instanceof NumberColumn && (mode === ERenderMode.CELL);
        return mode === ERenderMode.CELL;
    };
    MyLineChartRenderer.prototype.create = function (col) {
        return {
            template: "<div class=\"svg-container\"><svg class=\"svg-content\" preserveAspectRatio=\"xMidYMid meet\"><g><path class=\"areaChart\"></path><path class=\"lineChart\" fill=\"none\" stroke=\"".concat(base_color, "\" stroke-width=\"0.02px\"></path><g><line class=\"focus-line\"></line></g><g><text style=\"font-size:0.2px;\" class=\"focus-text\"></text></g><rect class=\"hover-rect\"></rect></g></svg></div>"),
            update: function (n, d) {
                if (renderMissingDOM(n, col, d)) {
                    return;
                }
                // get data
                var row = col.getMap(d);
                var data_mean_list = row[0]["value"];
                var data_var_list = row[1]["value"];
                // const data_max = col.getMax();
                var data_max = d3v5.max(data_mean_list, function (d) {
                    return +d;
                }) + 1;
                // const data_min = col.getMin();
                var data_min = d3v5.min(data_mean_list, function (d) {
                    return +d;
                }) - 1;
                // this is the ratio that the chart should have
                var rel_width = WIDTH_HEIGHT_RATIO; //data_mean_list.length/4;
                var rel_height = 1;
                var div = d3v5.select(n);
                var svg = div.select("svg");
                svg.attr("viewBox", "0 0 ".concat(rel_width, " ").concat(rel_height));
                // define x and y scales
                var x = d3v5
                    .scaleTime()
                    .domain([0, data_mean_list.length])
                    .range([0, rel_width]);
                var y = d3v5
                    .scaleLinear()
                    .domain([data_min, data_max])
                    .range([rel_height, 0]);
                // Show confidence interval
                svg
                    .select(".areaChart")
                    .datum(data_var_list)
                    .attr("fill", "#c1c1c14d")
                    .attr("stroke", "none")
                    .attr("d", d3v5
                    .area()
                    .x(function (d, i) {
                    return x(i);
                })
                    .y0(function (d, i) {
                    return y(data_mean_list[i] - d);
                })
                    .y1(function (d, i) {
                    return y(data_mean_list[i] + d);
                }));
                // draw the line chart
                var path = svg.select(".lineChart");
                path.datum(data_mean_list).attr("d", d3v5
                    .line()
                    .x(function (d, i) {
                    return x(i);
                }) // i/data_list.length
                    .y(function (d) {
                    return y(d);
                }) // 1-(d/data_max)
                );
                // add tooltips
                // https://www.d3-graph-gallery.com/graph/line_cursor.html
                // This allows to find the closest X index of the mouse:
                // var bisect = d3v5.bisector(function(d, i) { return i; }).left;
                // Create the line that travels along the x-axis of chart
                var focus = svg
                    .select(".focus-line")
                    .style("fill", "none")
                    .attr("stroke", "black")
                    .attr("stroke-width", "1%")
                    .attr("y1", "0")
                    .attr("y2", rel_height)
                    .attr("x1", "0")
                    .attr("x2", "0")
                    .style("opacity", 0);
                // Create the text that travels along the curve of chart
                var focusText = svg
                    .select(".focus-text")
                    .style("opacity", 0)
                    .attr("text-anchor", "left")
                    .attr("alignment-baseline", "middle")
                    .attr("letter-spacing", "0px");
                // Create a rect on top of the svg area: this rectangle recovers mouse position
                svg
                    .select(".hover-rect")
                    .style("fill", "none")
                    .style("pointer-events", "all")
                    .attr("width", "100%")
                    .attr("height", "100%")
                    .on("mouseover", mouseover)
                    .on("mousemove", mousemove)
                    .on("mouseout", mouseout);
                // What happens when the mouse move -> show the annotations at the right positions.
                function mouseover() {
                    focus.style("opacity", 1);
                    focusText.style("opacity", 1);
                }
                function mousemove() {
                    // recover coordinate we need
                    //@ts-ignore
                    var x0 = d3v5.mouse(this)[0];
                    // var y0 = d3v5.mouse(this)[1];
                    // var i = bisect(data, x0, 1);
                    // var x0 = x.invert(d3v5.mouse(this)[0]);
                    var i = Math.round(x.invert(x0)); // x0*data_list.length
                    i = Math.max(i, 0);
                    i = Math.min(data_mean_list.length - 1, i);
                    focus
                        .attr("x1", x(i)) // i/data_list.length
                        .attr("x2", x(i)); // i/data_list.length
                    // // position the text in a way that it is always readable
                    // if(x0 > rel_width/2){
                    //     x0 = x0-rel_width/2;
                    // }
                    focusText
                        .html("<tspan x='0' dy='1.2em'>step: " +
                        i +
                        "</tspan><tspan x='0' dy='1.2em'>mean: " +
                        Math.round(data_mean_list[i] * 100) / 100 +
                        "</tspan><tspan x='0' dy='1.2em'>var: " +
                        Math.round(data_var_list[i] * 100) / 100 +
                        "</tspan>")
                        .attr("x", 0) //x0
                        .attr("y", 0); //y0
                }
                function mouseout() {
                    focus.style("opacity", 0);
                    focusText.style("opacity", 0);
                }
            },
        };
    };
    return MyLineChartRenderer;
}());
export { MyLineChartRenderer };
var MySmilesStructureRenderer = /** @class */ (function () {
    function MySmilesStructureRenderer() {
        this.title = "Compound Structure";
        // better with background image, because with img tag the user might drag the img when they actually want to select several rows
        this.template = '<div style="background-size: contain; background-position: center; background-repeat: no-repeat;"></div>';
    }
    MySmilesStructureRenderer.prototype.canRender = function (col, mode) {
        return (col instanceof StructureImageColumn &&
            (mode === ERenderMode.CELL || mode === ERenderMode.GROUP));
    };
    MySmilesStructureRenderer.prototype.create = function (col) {
        return {
            template: this.template,
            update: function (n, d) {
                // @ts-ignore
                var smiles = d.v[col.desc.column];
                // TODO: Change undefined to id, how to get dataset id?
                CIMEBackendFromEnv.getStructureFromSmiles(undefined, smiles, false, null).then(function (x) {
                    if (x.length > 100) {
                        // check if it is actually long enogh to be an img
                        n.style.backgroundImage = "url('data:image/jpg;base64,".concat(x, "')");
                    }
                    else {
                        n.innerHTML = x;
                    }
                    n.alt = smiles;
                });
            },
        };
    };
    MySmilesStructureRenderer.prototype.createGroup = function (col, context) {
        return {
            template: this.template,
            update: function (n, group) {
                var formData = new FormData();
                return context.tasks
                    .groupRows(col, group, "string", function (rows) {
                    rows.every(function (row) {
                        var v = col.getLabel(row);
                        formData.append("smiles_list", v);
                        return true;
                    });
                })
                    .then(function () {
                    CIMEBackendFromEnv.getMCSFromSmilesList(formData).then(function (x) {
                        n.style.backgroundImage = "url('data:image/jpg;base64,".concat(x, "')");
                        n.alt = formData.getAll("smiles_list").toString();
                    });
                });
            },
        };
    };
    return MySmilesStructureRenderer;
}());
export { MySmilesStructureRenderer };
//# sourceMappingURL=LineUpContext.js.map