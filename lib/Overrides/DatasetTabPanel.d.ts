/// <reference types="react" />
import { ConnectedProps } from "react-redux";
export declare const LoadingIndicatorView: (props: any) => JSX.Element;
export declare const LoadingIndicatorDialog: (props: any) => JSX.Element;
declare const connector: import("react-redux").InferableComponentEnhancerWithProps<{
    setDataset: (value: any) => any;
}, {}>;
type Props = ConnectedProps<typeof connector>;
export declare const DatasetTabPanel: import("react-redux").ConnectedComponent<({ setDataset }: Props) => JSX.Element, import("react-redux").Omit<{
    setDataset: (value: any) => any;
}, "setDataset">>;
export {};
//# sourceMappingURL=DatasetTabPanel.d.ts.map