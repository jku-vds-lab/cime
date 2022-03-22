/// <reference types="react" />
import { ConnectedProps } from "react-redux";
import "./LineUpContext.scss";
import { IStringFilter, ERenderMode, IRenderContext, ICellRenderer, ICellRendererFactory, IDataRow, IGroupCellRenderer, StringColumn } from "lineupjs";
import { TestColumn } from "./LineUpClasses/TestColumn";
/**
 * Factory method which is declared here so we can get a static type in 'ConnectedProps'
 */
declare const connector: import("react-redux").InferableComponentEnhancerWithProps<{
    dataset: import("projection-space-explorer").Dataset;
    lineUpInput: import("../State/LineUpInputDuck").LineUpType;
    lineUpInput_data: import("projection-space-explorer").IVector[];
    lineUpInput_columns: {
        [name: string]: {
            distinct: any;
            isNumeric: boolean;
            metaInformation: any;
            featureType: import("projection-space-explorer").FeatureType;
            range: any;
            featureLabel: string;
            project: boolean;
        };
    };
    currentAggregation: {
        aggregation: number[];
        selectedClusters: import("@reduxjs/toolkit").EntityId[];
        source: "sample" | "cluster";
    };
    activeStory: import("projection-space-explorer").IBook;
    pointColorScale: any;
    channelColor: any;
    detailView: {
        open: boolean;
        active: string;
    };
} & {
    setCurrentAggregation: (samples: number[]) => any;
    setLineUpInput_visibility: (visibility: any) => any;
    setLineUpInput_lineup: (input: any) => any;
    setHoverstate: (state: any, updater: any) => any;
}, {}>;
/**
 * Type that holds the props we declared above in mapStateToProps and mapDispatchToProps
 */
declare type PropsFromRedux = ConnectedProps<typeof connector>;
/**
 * Type that holds every property that is relevant to our component, that is the props declared above + our OWN component props
 */
declare type Props = PropsFromRedux & {
    onFilter: any;
};
/**
 * Our component definition, by declaring our props with 'Props' we have static types for each of our property
 */
export declare const LineUpContext: import("react-redux").ConnectedComponent<({ dataset, lineUpInput, lineUpInput_data, lineUpInput_columns, currentAggregation, channelColor, setCurrentAggregation, setLineUpInput_lineup, setLineUpInput_visibility, onFilter, activeStory, pointColorScale, setHoverstate, detailView, }: Props) => JSX.Element, import("react-redux").Omit<{
    dataset: import("projection-space-explorer").Dataset;
    lineUpInput: import("../State/LineUpInputDuck").LineUpType;
    lineUpInput_data: import("projection-space-explorer").IVector[];
    lineUpInput_columns: {
        [name: string]: {
            distinct: any;
            isNumeric: boolean;
            metaInformation: any;
            featureType: import("projection-space-explorer").FeatureType;
            range: any;
            featureLabel: string;
            project: boolean;
        };
    };
    currentAggregation: {
        aggregation: number[];
        selectedClusters: import("@reduxjs/toolkit").EntityId[];
        source: "sample" | "cluster";
    };
    activeStory: import("projection-space-explorer").IBook;
    pointColorScale: any;
    channelColor: any;
    detailView: {
        open: boolean;
        active: string;
    };
} & {
    setCurrentAggregation: (samples: number[]) => any;
    setLineUpInput_visibility: (visibility: any) => any;
    setLineUpInput_lineup: (input: any) => any;
    setHoverstate: (state: any, updater: any) => any;
} & {
    onFilter: any;
}, "lineUpInput" | "dataset" | "setCurrentAggregation" | "setHoverstate" | "lineUpInput_data" | "lineUpInput_columns" | "currentAggregation" | "channelColor" | "setLineUpInput_lineup" | "setLineUpInput_visibility" | "activeStory" | "pointColorScale" | "detailView">>;
export interface IStructureFilter extends IStringFilter {
    filter: string;
    valid: Set<string>;
}
export declare class StructureImageColumn extends StringColumn {
    protected structureFilter: IStructureFilter | null;
    filter(row: IDataRow): boolean;
    isFiltered(): boolean;
    getFilter(): IStructureFilter;
    setFilter(filter: IStructureFilter | null): void;
}
export declare class MyLineChartRenderer implements ICellRendererFactory {
    readonly title: string;
    canRender(col: TestColumn, mode: ERenderMode): boolean;
    create(col: TestColumn): ICellRenderer;
}
export declare class MySmilesStructureRenderer implements ICellRendererFactory {
    readonly title: string;
    readonly template = "<div style=\"background-size: contain; background-position: center; background-repeat: no-repeat;\"></div>";
    canRender(col: StructureImageColumn, mode: ERenderMode): boolean;
    create(col: StructureImageColumn): ICellRenderer;
    createGroup(col: StructureImageColumn, context: IRenderContext): IGroupCellRenderer;
}
export {};
//# sourceMappingURL=LineUpContext.d.ts.map